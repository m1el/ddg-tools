//! Differential test against the Lean oracle (`lean/oracle`, command `rc`), SPEC §13.6.
//!
//! The oracle is an independent Float implementation, written from `SPEC.md` §2-§6 without this
//! crate's code, of: the energy and its differential (tasks `energy`), the dense Gram matrices `B`,
//! `B⁰`, `A` by their increments and the kernel matrices `K^B`, `K⁰`, at both quadratures (`gram`),
//! the constraint values and Jacobians of every constraint type together with one projected
//! Sobolev gradient from the dense saddle system (`saddle`), and one chord projection of a trial
//! step `γ − τ g̃` (`project`). This test is written from the oracle's README section "Repulsive
//! curves" and this crate's public API.
//!
//! * `fixtures_are_reproduced` (always runs): replays the committed golden fixtures
//!   `lean/oracle/fixtures/rc/*.in|out`.
//! * `fuzz_against_the_oracle` (`#[ignore]`): asks the built oracle for fresh seeded cases
//!   (`ddg-oracle rc-gen SEED COUNT FIRST`), has the oracle answer them (`ddg-oracle rc`) and
//!   replays every answer, in parallel over all cores. Run
//!
//!   ```text
//!   (cd lean/oracle && lake build)
//!   DDG_RC_FUZZ_CASES=20000 DDG_RC_FUZZ_SEED=1 \
//!     cargo test --release -p ddg-repulsive-curves --test oracle -- --ignored --nocapture
//!   ```
//!
//!   Environment: `DDG_RC_FUZZ_CASES` (default 4000), `DDG_RC_FUZZ_SEED` (default 1),
//!   `DDG_RC_FUZZ_THREADS` (default: all cores), `DDG_ORACLE` (path of the executable, default
//!   `lean/oracle/.lake/build/bin/ddg-oracle`), `DDG_RC_FUZZ_DIR` (where disagreeing inputs are
//!   written, default the target tmp dir).
//!
//! # Tolerances
//!
//! Both sides evaluate the same formulas in `f64` with different summation orders; tolerances are
//! normwise (relative to the largest entry of the object compared) because the entries of a
//! gradient or a Gram matrix span many orders of magnitude and the small ones carry the rounding
//! of the large ones. The observed worst ratios are printed by every run (`worst ratio`).
//!
//! | quantity | tolerance |
//! |---|---|
//! | energy | relative `1e-12` |
//! | differential | normwise `1e-11` (+ `1e-14` of the largest entry as floor, entries along the held axis `1e-12`) |
//! | `B`, `B⁰`, `A`, `K^B`, `K⁰` | normwise `1e-12` |
//! | constraint values `Φ` | `1e-12 (ℓ̄ + |Φ|)` |
//! | Jacobian `C` | normwise `1e-12` |
//! | `g̃` | normwise `1e-10` plus `1e-12 ℓ̄` (it vanishes at symmetric critical points) |
//! | `λ`, `dF·g̃`, `‖g̃‖` | normwise `1e-9` plus `1e-12` of the scale |
//! | all solve quantities (`g̃`, `λ`, `dF·g̃`, projection) | multiplied by `max(1, 50 eps κ / 1e-10)`, `κ` the condition number printed by the oracle (`cond`), for the few percent-of-a-percent ill-conditioned saddles |
//! | projection: status, number of iterations | exactly, unless a violation history entry is within `1e-5` of 1 (a decision at the tolerance) or the status is a rounding race |
//! | projection: history `max|Φ_r|/tol_r` | `1e-6` relative plus `1e-3` absolute (entries below 1 are rounding-level residuals) |
//! | projection: final positions | `1e-9 ℓ̄` |
//! | invalid input | the oracle names the error variant; this crate must return an error of the same variant (singular saddle: any of `SingularMetric`, `SingularSaddle`, `RedundantConstraints`) |
//!
//! Not modelled by the oracle (so never generated): intersecting edges (`EdgesIntersect`), tangent
//! constraints in planar mode (the oracle answers `Unsupported`, which is skipped), row dropping
//! of D9 (the oracle reports a singular saddle where this crate drops rows).

#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Arc;

use ddg_repulsive_curves::constraints::{
    Constraint, ConstraintSet, ImplicitSurface, TangentSite,
};
use ddg_repulsive_curves::energy::exact::make_pool;
use ddg_repulsive_curves::metric::sobolev_dense::GramGeometry;
use ddg_repulsive_curves::metric::{GramQuadrature, Metric, gram, kernel_matrices};
use ddg_repulsive_curves::solve::{
    ProjectionOptions, ProjectionStatus, Tolerance, dense_sobolev_solver, project_constraints,
    projected_gradient,
};
use ddg_repulsive_curves::{
    Axis, CurveNetwork, EnergyEval, Exponents, RcError, check_planar, energy,
};

type P = [f64; 3];

// ----- tolerances (see the table in the module docs) -----

const ENERGY_REL: f64 = 1e-12;
const DIFF_REL: f64 = 1e-11;
const DIFF_FLOOR: f64 = 1e-14;
const GRAM_REL: f64 = 1e-12;
const PHI_REL: f64 = 1e-12;
const JAC_REL: f64 = 1e-12;
const G_REL: f64 = 1e-10;
const G_ABS_LBAR: f64 = 1e-12;
const SCALAR_REL: f64 = 1e-9;
const HIST_REL: f64 = 1e-6;
const HIST_ABS: f64 = 1e-3;
const POS_LBAR: f64 = 1e-9;
const DECISION: f64 = 1e-5;
/// A saddle whose smallest/largest LU pivot is below this is at the edge of singular: if the oracle
/// calls it singular and this crate does not (or the reverse), neither is wrong.
const SINGULAR_AMBIGUITY: f64 = 1e-11;
/// The oracle prints the infinity-norm condition number `κ` of its equilibrated saddle matrix. A
/// solve is trusted to about `50 eps κ` (`eps = 1.1e-16`); where that exceeds the base tolerance of
/// the solve quantities (`1e-10`), the tolerances are multiplied by `50 eps κ / 1e-10`.
const COND_ERROR: f64 = 50.0 * 1.1e-16;
const SOLVE_BASE: f64 = 1e-10;

// ----- the text format -----

type Line = Vec<String>;

fn float(token: &str) -> f64 {
    assert_eq!(token.len(), 16, "a float is 16 hexadecimal digits: {token}");
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|e| panic!("{token}: {e}")))
}

fn int(token: &str) -> usize {
    token.parse().unwrap_or_else(|e| panic!("{token}: {e}"))
}

/// The `case NAME ... end` blocks of a stream.
fn blocks(text: &str) -> Vec<(String, Vec<Line>)> {
    let mut out = Vec::new();
    let mut current: Option<(String, Vec<Line>)> = None;
    for raw in text.lines() {
        let line: Line = raw.split_whitespace().map(str::to_string).collect();
        if line.is_empty() || line[0].starts_with('#') {
            continue;
        }
        match line[0].as_str() {
            "case" => current = Some((line[1].clone(), Vec::new())),
            "end" => out.extend(current.take()),
            _ => {
                if let Some((_, lines)) = current.as_mut() {
                    lines.push(line)
                }
            }
        }
    }
    out
}

#[derive(Debug, Clone)]
enum Spec {
    Barycenter(Option<P>),
    TotalLength(Option<f64>),
    EdgeLengths(Option<Vec<u32>>),
    Point(u32, Option<P>),
    Sphere(P, f64, Option<Vec<u32>>),
    Tangent(bool, u32, Option<P>), // (is_vertex, index, direction)
}

struct Sphere {
    c: P,
    r: f64,
}

impl ImplicitSurface for Sphere {
    /// `|x − c|² − r²`, the oracle's surface.
    fn value(&self, x: P) -> f64 {
        let d = [x[0] - self.c[0], x[1] - self.c[1], x[2] - self.c[2]];
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - self.r * self.r
    }
    fn gradient(&self, x: P) -> P {
        [
            2.0 * (x[0] - self.c[0]),
            2.0 * (x[1] - self.c[1]),
            2.0 * (x[2] - self.c[2]),
        ]
    }
}

struct Input {
    name: String,
    task: String,
    alpha: f64,
    beta: f64,
    quad: GramQuadrature,
    planar: Option<usize>,
    positions: Vec<P>,
    edges: Vec<[u32; 2]>,
    cons: Vec<Spec>,
    tol_mean: bool,
    tol_val: f64,
    maxit: u32,
}

fn p3(l: &Line, i: usize) -> P {
    [float(&l[i]), float(&l[i + 1]), float(&l[i + 2])]
}

fn ids(l: &Line, from: usize) -> Vec<u32> {
    l[from..].iter().map(|t| int(t) as u32).collect()
}

fn parse_input(name: &str, lines: &[Line]) -> Input {
    let mut input = Input {
        name: name.to_string(),
        task: String::new(),
        alpha: 3.0,
        beta: 6.0,
        quad: GramQuadrature::Midpoint,
        planar: None,
        positions: Vec::new(),
        edges: Vec::new(),
        cons: Vec::new(),
        tol_mean: true,
        tol_val: 1e-4,
        maxit: 3,
    };
    for l in lines {
        match l[0].as_str() {
            "task" => input.task = l[1].clone(),
            "exponents" => {
                input.alpha = float(&l[1]);
                input.beta = float(&l[2]);
            }
            "quadrature" => {
                input.quad = if l[1] == "trapezoid" {
                    GramQuadrature::Trapezoid
                } else {
                    GramQuadrature::Midpoint
                }
            }
            "planar" => input.planar = Some(int(&l[1])),
            "matrices" | "step" | "steprel" => {}
            "tol" => {
                input.tol_mean = l[1] == "mean";
                input.tol_val = float(&l[2]);
            }
            "maxit" => input.maxit = int(&l[1]) as u32,
            "p" => input.positions.push(p3(l, 1)),
            "e" => input.edges.push([int(&l[1]) as u32, int(&l[2]) as u32]),
            "net" => {}
            "c" => input.cons.push(match l[1].as_str() {
                "barycenter" => Spec::Barycenter((l.len() == 5).then(|| p3(l, 2))),
                "total-length" => Spec::TotalLength((l.len() == 3).then(|| float(&l[2]))),
                "edge-lengths" => {
                    Spec::EdgeLengths((l.len() < 3 || l[2] != "all").then(|| ids(l, 2)))
                }
                "point" => Spec::Point(int(&l[2]) as u32, (l.len() == 6).then(|| p3(l, 3))),
                "sphere" => Spec::Sphere(
                    p3(l, 2),
                    float(&l[5]),
                    (l[6] != "all").then(|| ids(l, 6)),
                ),
                "tangent" => Spec::Tangent(
                    l[2] == "vertex",
                    int(&l[3]) as u32,
                    (l.len() == 7).then(|| p3(l, 4)),
                ),
                other => panic!("unknown constraint {other}"),
            }),
            other => panic!("unknown input line {other}"),
        }
    }
    input
}

fn constraints_of(input: &Input) -> Vec<Constraint> {
    input
        .cons
        .iter()
        .map(|s| match s {
            Spec::Barycenter(t) => Constraint::Barycenter { target: *t },
            Spec::TotalLength(t) => Constraint::TotalLength { target: *t },
            Spec::EdgeLengths(e) => Constraint::EdgeLengths {
                edges: e.clone(),
                targets: None,
            },
            Spec::Point(v, t) => Constraint::Point {
                vertex: *v,
                target: *t,
            },
            Spec::Sphere(c, r, vs) => Constraint::Surface {
                vertices: vs.clone(),
                surface: Arc::new(Sphere { c: *c, r: *r }),
            },
            Spec::Tangent(vertex, i, d) => Constraint::Tangent {
                site: if *vertex {
                    TangentSite::Vertex(*i)
                } else {
                    TangentSite::Edge(*i)
                },
                direction: *d,
            },
        })
        .collect()
}

// ----- statistics -----

#[derive(Default)]
struct Stats {
    cases: usize,
    by_task: BTreeMap<String, usize>,
    ok_cases: usize,
    error_cases: usize,
    skipped: usize,
    tie_skips: usize,
    ill_conditioned: usize,
    ambiguous_singular: usize,
    statuses: BTreeMap<String, usize>,
    comparisons: usize,
    /// per quantity: (largest ratio error/tolerance, number of comparisons)
    worst: BTreeMap<&'static str, (f64, usize)>,
    problems: Vec<(String, String)>,
}

impl Stats {
    fn merge(&mut self, o: Stats) {
        self.cases += o.cases;
        for (k, v) in o.by_task {
            *self.by_task.entry(k).or_default() += v;
        }
        self.ok_cases += o.ok_cases;
        self.error_cases += o.error_cases;
        self.skipped += o.skipped;
        self.tie_skips += o.tie_skips;
        self.ill_conditioned += o.ill_conditioned;
        self.ambiguous_singular += o.ambiguous_singular;
        for (k, v) in o.statuses {
            *self.statuses.entry(k).or_default() += v;
        }
        self.comparisons += o.comparisons;
        for (k, (r, n)) in o.worst {
            let e = self.worst.entry(k).or_insert((0.0, 0));
            e.0 = e.0.max(r);
            e.1 += n;
        }
        self.problems.extend(o.problems);
    }

    fn problem(&mut self, name: &str, msg: String) {
        self.problems.push((name.to_string(), msg));
    }

    /// Records `ratio = error / tolerance`; a problem if it exceeds 1 or is NaN.
    fn check(&mut self, name: &str, what: &'static str, ratio: f64, detail: impl FnOnce() -> String) {
        self.comparisons += 1;
        let e = self.worst.entry(what).or_insert((0.0, 0));
        e.1 += 1;
        if ratio.is_nan() || ratio > e.0 {
            e.0 = if ratio.is_nan() { f64::INFINITY } else { ratio };
        }
        if !(ratio <= 1.0) {
            let d = detail();
            self.problem(name, format!("{what}: ratio {ratio:.3e} (> 1): {d}"));
        }
    }

    fn report(&self) -> String {
        let mut s = format!(
            "{} cases ({} answered, {} errors, {} skipped, {} decision ties, {} ill-conditioned saddles with scaled tolerances, {} singular-or-not ambiguities), {} comparisons, {} disagreements\n",
            self.cases,
            self.ok_cases,
            self.error_cases,
            self.skipped,
            self.tie_skips,
            self.ill_conditioned,
            self.ambiguous_singular,
            self.comparisons,
            self.problems.len()
        );
        s += &format!("  by task: {:?}\n", self.by_task);
        s += &format!("  projection endings: {:?}\n", self.statuses);
        for (k, (r, n)) in &self.worst {
            s += &format!("  worst ratio error/tolerance {k:>26}: {r:.3e}  ({n} comparisons)\n");
        }
        s
    }
}

// ----- comparison helpers -----

fn max_abs(a: &[f64]) -> f64 {
    a.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// `max|a - b| / (rel * max(max|a|, max|b|) + floor)`; NaN-safe (NaN gives NaN).
fn normwise(a: &[f64], b: &[f64], rel: f64, floor: f64) -> f64 {
    if a.len() != b.len() {
        return f64::INFINITY;
    }
    let mut err = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        let d = (x - y).abs();
        if d.is_nan() {
            return f64::NAN;
        }
        err = err.max(d);
    }
    err / (rel * max_abs(a).max(max_abs(b)) + floor)
}

fn flat(v: &[P]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn rows_with<'a>(out: &'a [Line], tag: &str) -> Vec<&'a Line> {
    out.iter().filter(|l| l[0] == tag).collect()
}

fn first<'a>(out: &'a [Line], tag: &str) -> Option<&'a Line> {
    out.iter().find(|l| l[0] == tag)
}

/// Oracle rows `tag i v...` as a flat row-major matrix.
fn oracle_matrix(rows: &[&Line]) -> Vec<f64> {
    rows.iter()
        .flat_map(|l| l[2..].iter().map(|t| float(t)))
        .collect()
}

fn oracle_points(out: &[Line], tag: &str) -> Vec<P> {
    rows_with(out, tag)
        .iter()
        .map(|l| p3(l, 2))
        .collect()
}

// ----- the Rust side -----

fn exponents(input: &Input) -> Result<Exponents, RcError> {
    Exponents::new(input.alpha, input.beta)
}

fn network(input: &Input) -> Result<CurveNetwork, RcError> {
    let net = CurveNetwork::new(input.positions.len(), input.edges.clone())?;
    if let Some(a) = input.planar {
        let axis = [Axis::X, Axis::Y, Axis::Z][a];
        check_planar(&input.positions, axis)?;
    }
    net.validate_positions(&input.positions)?;
    Ok(net)
}

fn axis_of(input: &Input) -> Option<Axis> {
    input.planar.map(|a| [Axis::X, Axis::Y, Axis::Z][a])
}

fn variant(e: &RcError) -> String {
    format!("{e:?}")
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .unwrap_or("")
        .to_string()
}

fn error_matches(oracle_kind: &str, e: &RcError) -> bool {
    let v = variant(e);
    match oracle_kind {
        "SingularSaddle" => [
            "SingularMetric",
            "SingularSaddle",
            "RedundantConstraints",
            "TranslationNotFixed",
        ]
        .contains(&v.as_str()),
        k => k == v,
    }
}

/// Everything the Rust side computes for the two linear-algebra tasks.
struct Solved {
    phi: Vec<f64>,
    c: Vec<f64>,
    n_rows: usize,
    n_dofs: usize,
    g: Vec<P>,
    lambda: Vec<f64>,
    dfg: f64,
    gnorm: f64,
    mean_edge: f64,
    diff_max: f64,
    layout_coords: Vec<usize>,
}

fn solve_rust(input: &Input, net: &CurveNetwork, ex: Exponents) -> Result<(Solved, ConstraintSet, ddg_repulsive_curves::solve::DenseSaddle), RcError> {
    let pos = &input.positions;
    let pool = Arc::new(make_pool(Some(1))?);
    let report = energy(net, pos, ex, &EnergyEval::Exact, Some(1))?;
    let cons = constraints_of(input);
    let cs = ConstraintSet::new(net, pos, &cons, axis_of(input))?;
    let phi = cs.values(net, pos)?;
    let jac = cs.jacobian(net, pos)?;
    let (k, nd) = (jac.nrows(), jac.ncols());
    let mut c = vec![0.0; k * nd];
    for (r, col, v) in jac.triplets() {
        c[r as usize * nd + col as usize] += v;
    }
    let mut solver = dense_sobolev_solver(pool, net, pos, ex, input.quad, &cs)?;
    let layout = cs.layout();
    let pg = projected_gradient(&mut solver, layout, &report.differential)?;
    let mean_edge = net.mean_edge_length(pos);
    Ok((
        Solved {
            phi,
            c,
            n_rows: k,
            n_dofs: nd,
            g: pg.direction,
            lambda: pg.lambda,
            dfg: pg.sobolev_norm_sq,
            gnorm: pg.gradient_norm,
            mean_edge,
            diff_max: max_abs(&flat(&report.differential)),
            layout_coords: layout.coords().to_vec(),
        },
        cs,
        solver,
    ))
}

// ----- replay -----

fn replay(input: &Input, out: &[Line], stats: &mut Stats) {
    stats.cases += 1;
    *stats.by_task.entry(input.task.clone()).or_default() += 1;
    let name = input.name.as_str();
    let status = first(out, "status").map(|l| l[1].as_str()).unwrap_or("?");
    if status == "error" {
        let l = first(out, "status").unwrap();
        let kind = l.get(2).cloned().unwrap_or_default();
        if kind == "Unsupported" {
            stats.skipped += 1;
            return;
        }
        stats.error_cases += 1;
        match rust_pipeline_error(input) {
            (None, Some(ratio)) if kind == "SingularSaddle" && ratio < SINGULAR_AMBIGUITY => {
                // the saddle is at the edge of singular for both: either answer is right
                stats.ambiguous_singular += 1;
            }
            (Some(e), _) => {
                if !error_matches(&kind, &e) {
                    stats.problem(
                        name,
                        format!("the oracle reports {kind} ({}); this crate returns {e:?}", l[3..].join(" ")),
                    );
                }
            }
            (None, _) => stats.problem(
                name,
                format!("the oracle reports {kind} ({}); this crate accepts the input", l[3..].join(" ")),
            ),
        }
        stats.comparisons += 1;
        return;
    }
    if status != "ok" {
        stats.problem(name, "no status line from the oracle".into());
        return;
    }
    stats.ok_cases += 1;
    let ex = match exponents(input) {
        Ok(e) => e,
        Err(e) => return stats.problem(name, format!("the oracle accepts the exponents; this crate: {e:?}")),
    };
    let net = match network(input) {
        Ok(n) => n,
        Err(e) => return stats.problem(name, format!("the oracle accepts the network; this crate: {e:?}")),
    };
    match input.task.as_str() {
        "energy" => replay_energy(input, &net, ex, out, stats),
        "gram" => replay_gram(input, &net, ex, out, stats),
        _ => replay_solve(input, &net, ex, out, stats),
    }
}

/// The error this crate returns for the input, running the same pipeline as the task.
fn rust_pipeline_error(input: &Input) -> (Option<RcError>, Option<f64>) {
    let ex = match exponents(input) {
        Ok(e) => e,
        Err(e) => return (Some(e), None),
    };
    let net = match network(input) {
        Ok(n) => n,
        Err(e) => return (Some(e), None),
    };
    match input.task.as_str() {
        "energy" | "gram" => (None, None),
        _ => match solve_rust(input, &net, ex) {
            Ok((_, _, solver)) => (None, Some(solver.pivot_ratio())),
            Err(e) => (Some(e), None),
        },
    }
}

fn replay_energy(input: &Input, net: &CurveNetwork, ex: Exponents, out: &[Line], stats: &mut Stats) {
    let name = input.name.as_str();
    let rep = match energy(net, &input.positions, ex, &EnergyEval::Exact, Some(1)) {
        Ok(r) => r,
        Err(e) => return stats.problem(name, format!("energy: {e:?}")),
    };
    let oe = float(&first(out, "energy").expect("energy line")[1]);
    let ratio = (rep.energy - oe).abs() / (ENERGY_REL * rep.energy.abs().max(oe.abs()));
    stats.check(name, "energy", ratio, || format!("{} vs oracle {}", rep.energy, oe));
    let od = oracle_points(out, "d");
    let mine = rep.differential;
    // components along the held axis are zero in the oracle; this crate leaves rounding
    let (mut a, mut b) = (flat(&mine), flat(&od));
    if let Some(ax) = input.planar {
        let scale = max_abs(&b).max(max_abs(&a));
        for v in 0..mine.len() {
            if (mine[v][ax]).abs() > 1e-12 * scale.max(1e-300) {
                stats.problem(name, format!("differential along the held axis at vertex {v}: {}", mine[v][ax]));
            }
            a[3 * v + ax] = 0.0;
            b[3 * v + ax] = 0.0;
        }
    }
    let ratio = normwise(&a, &b, DIFF_REL, DIFF_FLOOR * max_abs(&b));
    stats.check(name, "differential", ratio, || {
        format!("max abs diff {:.3e} of {:.3e}", a.iter().zip(&b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs())), max_abs(&b))
    });
}

fn replay_gram(input: &Input, net: &CurveNetwork, ex: Exponents, out: &[Line], stats: &mut Stats) {
    let name = input.name.as_str();
    let n = net.n_vertices();
    let m = net.n_edges();
    let pool = make_pool(Some(1)).unwrap();
    let geom = GramGeometry::new(net, &input.positions, ex);
    // the sections in the oracle's order
    let mut sections: Vec<(&str, GramQuadrature, Vec<&Line>)> = Vec::new();
    let mut cur: Option<usize> = None;
    for l in out {
        if l[0] == "section" {
            let q = if l[1] == "midpoint" {
                GramQuadrature::Midpoint
            } else {
                GramQuadrature::Trapezoid
            };
            sections.push((if l[1] == "midpoint" { "midpoint" } else { "trapezoid" }, q, Vec::new()));
            cur = Some(sections.len() - 1);
        } else if let Some(c) = cur {
            sections[c].2.push(l);
        }
    }
    if sections.len() != 2 {
        return stats.problem(name, format!("expected two sections, found {}", sections.len()));
    }
    for (qname, q, lines) in sections {
        let g = match gram(net, &input.positions, ex, &Metric::Sobolev { quadrature: q }, Some(1)) {
            Ok(g) => g,
            Err(e) => return stats.problem(name, format!("gram: {e:?}")),
        };
        let pick = |tag: &str| -> Vec<&Line> { lines.iter().copied().filter(|l| l[0] == tag).collect() };
        let a = oracle_matrix(&pick("A"));
        stats.check(name, "gram A", normwise(g.a.as_slice(), &a, GRAM_REL, 0.0), || qname.into());
        let b = pick("B");
        if !b.is_empty() {
            stats.check(name, "gram B", normwise(g.b.as_slice(), &oracle_matrix(&b), GRAM_REL, 0.0), || qname.into());
            stats.check(name, "gram B0", normwise(g.b0.as_slice(), &oracle_matrix(&pick("B0")), GRAM_REL, 0.0), || qname.into());
            let km = kernel_matrices(&pool, net, &geom, q);
            let okb = oracle_matrix(&pick("KB"));
            let ok0 = oracle_matrix(&pick("K0"));
            assert_eq!(okb.len(), m * m);
            stats.check(name, "kernel K^B", normwise(km.kb.as_slice(), &okb, GRAM_REL, 0.0), || qname.into());
            stats.check(name, "kernel K0", normwise(km.k0.as_slice(), &ok0, GRAM_REL, 0.0), || qname.into());
        }
        let _ = n;
    }
}

fn replay_solve(input: &Input, net: &CurveNetwork, ex: Exponents, out: &[Line], stats: &mut Stats) {
    let name = input.name.as_str();
    let (s, cs, mut solver) = match solve_rust(input, net, ex) {
        Ok(x) => x,
        Err(e) => {
            let piv = float(&first(out, "pivot").unwrap()[1]);
            if piv < SINGULAR_AMBIGUITY && variant(&e).starts_with("Singular") {
                stats.ambiguous_singular += 1;
                return;
            }
            return stats.problem(name, format!("the oracle solves the saddle; this crate: {e:?}"));
        }
    };
    let cond = float(&first(out, "cond").unwrap()[1]);
    let mult = (COND_ERROR * cond / SOLVE_BASE).max(1.0);
    if mult > 1.0 {
        stats.ill_conditioned += 1;
    }
    let rows = first(out, "rows").expect("rows line");
    let (ok, ond) = (int(&rows[1]), int(&rows[2]));
    if ok != s.n_rows || ond != s.n_dofs {
        return stats.problem(name, format!("the saddle has {ok} rows and {ond} unknowns in the oracle, {} and {} here", s.n_rows, s.n_dofs));
    }
    let lbar = s.mean_edge;
    let phi_o: Vec<f64> = first(out, "phi").unwrap()[1..].iter().map(|t| float(t)).collect();
    let phi_err = s.phi.iter().zip(&phi_o).fold(0.0f64, |m, (a, b)| m.max((a - b).abs() / (PHI_REL * (lbar + a.abs().max(b.abs())))));
    if s.phi.len() != phi_o.len() {
        return stats.problem(name, "row count of Φ differs".into());
    }
    stats.check(name, "constraint values", phi_err, || format!("{:?} vs {:?}", s.phi, phi_o));
    if input.task == "saddle" {
        let co = oracle_matrix(&rows_with(out, "c"));
        stats.check(name, "Jacobian C", normwise(&s.c, &co, JAC_REL, 0.0), || "C".into());
        let g_o = oracle_points(out, "g");
        let g_mine = flat(&s.g);
        let g_oracle = flat(&g_o);
        stats.check(name, "projected gradient", normwise(&g_mine, &g_oracle, G_REL * mult, G_ABS_LBAR * lbar * mult), || {
            format!("max|g| {:.3e}", max_abs(&g_oracle))
        });
        let lam_o: Vec<f64> = first(out, "lambda").unwrap()[1..].iter().map(|t| float(t)).collect();
        let dfg_o = float(&first(out, "dfg").unwrap()[1]);
        let gn_o = float(&first(out, "gnorm").unwrap()[1]);
        // scales of the quantities: lambda by its own max, dF.g by the product of the norms
        let gscale = max_abs(&g_oracle).max(1e-300);
        stats.check(name, "multipliers", normwise(&s.lambda, &lam_o, SCALAR_REL * mult, 0.0), || "λ".into());
        let dscale = (dfg_o.abs()).max(1e-300);
        stats.check(name, "dF.g", (s.dfg - dfg_o).abs() / (SCALAR_REL * mult * dscale + s.n_dofs as f64 * s.diff_max * G_ABS_LBAR * lbar * mult), || format!("{} vs {}", s.dfg, dfg_o));
        stats.check(name, "gradient norm", (s.gnorm - gn_o).abs() / (SCALAR_REL * mult * gn_o.abs() + G_ABS_LBAR * lbar * mult), || format!("{} vs {}", s.gnorm, gn_o));
        let _ = gscale;
        return;
    }
    // project: build the trial from the oracle's g̃ and tau, so that the projection is tested alone
    let g_o = oracle_points(out, "g");
    let tau = float(&first(out, "tau").unwrap()[1]);
    let lbar_o = float(&first(out, "lbar").unwrap()[1]);
    stats.check(name, "mean edge length", (lbar - lbar_o).abs() / (1e-14 * lbar), || format!("{lbar} vs {lbar_o}"));
    let coords = s.layout_coords.clone();
    let trial: Vec<P> = input
        .positions
        .iter()
        .zip(&g_o)
        .map(|(p, g)| {
            let mut q = *p;
            for &c in &coords {
                q[c] -= tau * g[c];
            }
            q
        })
        .collect();
    let tolerance = if input.tol_mean {
        Tolerance::MeanEdgeLength(input.tol_val)
    } else {
        Tolerance::Absolute(input.tol_val)
    };
    let opts = ProjectionOptions {
        tolerance,
        max_iterations: input.maxit,
    };
    let pr = match project_constraints(&mut solver, &cs, net, &trial, &opts, lbar) {
        Ok(p) => p,
        Err(e) => return stats.problem(name, format!("projection: {e:?}")),
    };
    let res = first(out, "result").expect("result line");
    let o_iters = int(res.last().unwrap());
    let o_status = res[1..res.len() - 1].join(" ");
    let my_status = match pr.status {
        ProjectionStatus::Converged => "converged".to_string(),
        ProjectionStatus::MaxIterations => "max-iterations".to_string(),
        ProjectionStatus::Diverged => "diverged".to_string(),
        ProjectionStatus::TangentFlipped { constraint } => format!("tangent-flipped {constraint}"),
        _ => "other".to_string(),
    };
    *stats.statuses.entry(o_status.split(' ').next().unwrap().to_string()).or_default() += 1;
    let hist_o: Vec<f64> = first(out, "hist").unwrap()[1..].iter().map(|t| float(t)).collect();
    let near_tie = hist_o.iter().chain(pr.history.iter()).any(|h| (h - 1.0).abs() < DECISION);
    let n_common = hist_o.len().min(pr.history.len());
    // history entries up to the common length
    for i in 0..n_common {
        let (a, b) = (pr.history[i], hist_o[i]);
        let ratio = (a - b).abs() / ((HIST_REL * a.abs().max(b.abs()) + HIST_ABS) * mult);
        stats.check(name, "projection history", ratio, || format!("entry {i}: {a} vs {b} (all: {:?} vs {:?})", pr.history, hist_o));
    }
    if near_tie {
        stats.tie_skips += 1;
        return;
    }
    stats.comparisons += 1;
    if my_status != o_status || pr.iterations as usize != o_iters {
        return stats.problem(
            name,
            format!(
                "projection ends '{my_status}' after {} iterations here, '{o_status}' after {o_iters} in the oracle (histories {:?} vs {:?})",
                pr.iterations, pr.history, hist_o
            ),
        );
    }
    let x_o = oracle_points(out, "x");
    let ratio = {
        let (a, b) = (flat(&pr.positions), flat(&x_o));
        let err = a.iter().zip(&b).fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
        err / (POS_LBAR * lbar * mult)
    };
    stats.check(name, "projected positions", ratio, || "x".into());
    let mv_o = float(&first(out, "maxviol").unwrap()[1]);
    // the final violation is rounding-level after convergence; compare only unconverged endings
    if my_status != "converged" {
        stats.check(name, "final violation", (pr.max_violation - mv_o).abs() / ((1e-6 * mv_o.abs() + 1e-12 * lbar) * mult), || format!("{} vs {}", pr.max_violation, mv_o));
    }
}

fn replay_text(inputs: &str, outputs: &str) -> Stats {
    let mut stats = Stats::default();
    let outs: HashMap<String, Vec<Line>> = blocks(outputs).into_iter().collect();
    for (name, lines) in blocks(inputs) {
        let input = parse_input(&name, &lines);
        match outs.get(&name) {
            Some(out) => replay(&input, out, &mut stats),
            None => {
                stats.cases += 1;
                stats.problem(&name, "no answer from the oracle".into())
            }
        }
    }
    stats
}

// ----- fixtures -----

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/rc")
}

#[test]
fn fixtures_are_reproduced() {
    let dir = fixtures_dir();
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|e| {
            let name = e.unwrap().file_name().into_string().unwrap();
            name.strip_suffix(".in").map(str::to_string)
        })
        .collect();
    names.sort();
    assert!(names.len() >= 40, "only {} fixtures in {}", names.len(), dir.display());
    let mut total = Stats::default();
    for name in &names {
        let inputs = std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap();
        let outputs = std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap();
        total.merge(replay_text(&inputs, &outputs));
    }
    println!("{}", total.report());
    for (name, msg) in &total.problems {
        println!("DISAGREEMENT {name}: {msg}");
    }
    assert!(total.problems.is_empty(), "{} disagreements with the golden fixtures", total.problems.len());
}

// ----- fuzzing -----

fn oracle_path() -> PathBuf {
    match std::env::var("DDG_ORACLE") {
        Ok(p) => PathBuf::from(p),
        Err(_) => PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/.lake/build/bin/ddg-oracle"),
    }
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

/// Generate `count` cases starting at index `first`, answer them, replay them.
fn fuzz_chunk(oracle: &Path, seed: usize, first: usize, count: usize) -> (Stats, String) {
    let generated = Command::new(oracle)
        .args(["rc-gen", &seed.to_string(), &count.to_string(), &first.to_string()])
        .output()
        .expect("run rc-gen");
    assert!(generated.status.success(), "rc-gen failed");
    let inputs = String::from_utf8(generated.stdout).unwrap();
    let mut child = Command::new(oracle)
        .arg("rc")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("run rc");
    let mut stdin = child.stdin.take().unwrap();
    let text = inputs.clone();
    let writer = std::thread::spawn(move || {
        stdin.write_all(text.as_bytes()).unwrap();
    });
    let answer = child.wait_with_output().expect("wait for rc");
    writer.join().unwrap();
    assert!(answer.status.success(), "rc failed");
    let outputs = String::from_utf8(answer.stdout).unwrap();
    (replay_text(&inputs, &outputs), inputs)
}

#[test]
#[ignore = "needs the built oracle (cd lean/oracle && lake build); see the module docs"]
fn fuzz_against_the_oracle() {
    let oracle = oracle_path();
    assert!(oracle.exists(), "{} is not built", oracle.display());
    let cases = env_usize("DDG_RC_FUZZ_CASES", 4000);
    let seed = env_usize("DDG_RC_FUZZ_SEED", 1);
    let threads = env_usize(
        "DDG_RC_FUZZ_THREADS",
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
    );
    let chunk = 20;
    let n_chunks = cases.div_ceil(chunk);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let total = std::sync::Mutex::new(Stats::default());
    let inputs_of_problems = std::sync::Mutex::new(HashMap::<String, String>::new());
    let started = std::time::Instant::now();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let j = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                    if j >= n_chunks {
                        break;
                    }
                    let first = j * chunk;
                    let count = chunk.min(cases - first);
                    let (stats, inputs) = fuzz_chunk(&oracle, seed, first, count);
                    if !stats.problems.is_empty() {
                        let mut map = inputs_of_problems.lock().unwrap();
                        for (name, _) in &stats.problems {
                            if let Some((_, lines)) = blocks(&inputs).into_iter().find(|(n, _)| n == name) {
                                let mut text = format!("case {name}\n");
                                for l in lines {
                                    text += &(l.join(" ") + "\n");
                                }
                                text += "end\n";
                                map.insert(name.clone(), text);
                            }
                        }
                    }
                    total.lock().unwrap().merge(stats);
                }
            });
        }
    });
    let total = total.into_inner().unwrap();
    println!("seed {seed}, {cases} cases, {threads} threads, {:.1} s", started.elapsed().as_secs_f64());
    println!("{}", total.report());
    let dir = std::env::var("DDG_RC_FUZZ_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_TARGET_TMPDIR")));
    let map = inputs_of_problems.into_inner().unwrap();
    for (name, msg) in total.problems.iter().take(40) {
        println!("DISAGREEMENT {name}: {msg}");
    }
    if !map.is_empty() {
        std::fs::create_dir_all(&dir).unwrap();
        for (name, text) in &map {
            std::fs::write(dir.join(format!("{name}.in")), text).unwrap();
        }
        println!("inputs of the disagreeing cases written to {}", dir.display());
    }
    assert!(total.problems.is_empty(), "{} disagreements", total.problems.len());
}
