//! Differential test against the Lean oracle (`lean/oracle`, command `rs`), SPEC §13.5.
//!
//! The oracle is an independent Float implementation, written from `SPEC.md` alone, of: the exact
//! tangent-point energy and its gradient (`energy`), `D_f` (`df`), the dense metric `L_τ`, `B`,
//! `B⁰`, `A`, `W`, `M` (`metric`), constraint values and Jacobians (`constraints`), the nested
//! saddle with Schur direction and Newton correction (`saddle`), the dense preconditioner
//! (`precond`), one dense projected-gradient step (`step`) and one remeshing operation
//! (`remesh-ops`). This test is written from the oracle's README section "Repulsive surfaces" and
//! this crate's public API.
//!
//! * `fixtures_are_reproduced` (always runs): replays the golden fixtures
//!   `lean/oracle/fixtures/rs/*.in|out` (tasks `energy`, `df`, `metric`, `constraints`, `saddle`,
//!   `precond`, `remesh-ops`, and the error fixtures) and asserts a minimum count per task.
//! * `fuzz_against_the_oracle` (`#[ignore]`): asks the built oracle for fresh seeded cases
//!   (`ddg-oracle rs-gen SEED COUNT FIRST`), has it answer them (`ddg-oracle rs`) and replays every
//!   answer in parallel. Run
//!
//!   ```text
//!   (cd lean/oracle && lake build)
//!   DDG_RS_FUZZ_CASES=2000 DDG_RS_FUZZ_SEED=1 \
//!     cargo test --release -p ddg-repulsive-surfaces --test oracle -- --ignored --nocapture
//!   ```
//!
//!   Environment: `DDG_RS_FUZZ_CASES` (default 2000), `DDG_RS_FUZZ_SEED` (1), `DDG_RS_FUZZ_THREADS`
//!   (all cores), `DDG_ORACLE` (path of the executable, default
//!   `lean/oracle/.lake/build/bin/ddg-oracle`), `DDG_RS_FUZZ_DIR` (where the inputs of disagreeing
//!   cases are written, default the target tmp dir).
//!
//! # Tolerances
//!
//! Both sides evaluate the same formulas in `f64` with different summation orders. Every
//! comparison is normwise: `max |got − want| / max(max |want|, floor)`, because the entries of a
//! gradient or a matrix span orders of magnitude and the small ones carry the rounding of the large
//! ones. The oracle itself reports a check tolerance of `1e-9`. The worst ratio `error / tolerance`
//! per quantity is printed by every run.
//!
//! | quantity | tolerance |
//! |---|---|
//! | energy value | `1e-12` |
//! | energy gradient | `1e-12` |
//! | `D_f` (`ū`, gradient) | `1e-12` |
//! | `L_τ`, `B`, `B⁰`, `A`, `W`, `W M⁻¹ W`, `M` | `1e-12` |
//! | constraint values `Φ` | `1e-12` of `max(1, targets, |Φ|)` |
//! | Jacobian `C` | `1e-12` |
//! | saddle: `g` | `1e-12`; `x`, `⟨g,x⟩`, `Φ_N`, `S`, `h` | `1e-10` (solves; observed worst ~1.6e-11) |
//! | precond: `ε`, `L_{2−s}` | `1e-12`; `L′`, `P⁻¹ r` | `1e-10` |
//! | remesh-ops | connectivity, parents, vertex map exactly; positions `1e-14` |
//! | step | `F₀`, trial energies, `F₁`, `F_final`, `‖x‖`, `⟨g,x⟩`, `τ`: `1e-10`; `ℓ̄₀` `1e-12`; positions, `Φ`: `1e-10`; outcomes, halvings, Newton rounds exactly |
//! | error cases | the oracle names the kind; this crate must return the matching `RsError` variant (table in `error_matches`) |
//!
//! # Skips and special readings (all justified)
//!
//! * `rs-error-singular` (the same pin twice, saddle): the oracle never drops rows and reports
//!   `SingularSaddle`; this crate's rank check (SPEC D10) drops the duplicate pin. The test asserts
//!   that the Rust constraint set succeeds and reports a dropped row.
//! * `rs-mobius-flip` (and any `InconsistentOrientation` of a flip): the Rust flip reorders the
//!   two faces locally and succeeds (SPEC §10); the test asserts that it succeeds and reports
//!   `reordered`.
//! * `constraints` cases that leave a component's translation free: the oracle's translation check
//!   belongs to the saddle, this crate checks it when the constraint set is built; the test adds a
//!   pin on every component and compares the rows of the oracle only.
//! * `step` (one dense projected-gradient step with the standard Armijo rule, WP5): replayed with
//!   [`Session`] (`Backend::Dense`, `EnergyEval::Exact`, per-component barycenters, no remeshing,
//!   `StartRule::{Absolute, MaxDisplacement, Reference}`, `max_halvings`, `NewtonSteps`). The
//!   oracle's objective is `ℰ̂` with the metric `w_TP A`; the crate's is `w_TP ℰ̂` with `w_TP A`
//!   (the reference's reading of the weight, D6), whose trajectory is the same with `τ` scaled by
//!   `1/w_TP`: the replay passes `τ_abs / w_TP` for the absolute start and compares `F / w_TP`,
//!   `⟨g, x⟩ / w_TP²`, `‖x‖ / w_TP`, `τ w_TP`. Cases with a Willmore term in the metric but no
//!   Willmore energy (`w_W > 0`; the session adds the term only with a Willmore penalty) are
//!   skipped. If a trial's Armijo outcome differs within a margin of `1e-9` (rounding), the case is
//!   recorded as a special reading and its remaining comparisons are skipped.
//! * Oracle answers `Unsupported` are skipped.

#![allow(clippy::needless_range_loop)]

#[path = "wp2_fixture/mod.rs"]
mod fixture;

use std::collections::BTreeMap;
use std::io::Write;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use ddg_repulsive_surfaces::constraints::{Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::{energy_and_differential, make_pool};
use ddg_repulsive_surfaces::geometry::FaceGeometry;
use ddg_repulsive_surfaces::metric::dense::{
    dense_metric, fractional_operator, to_dense, willmore_matrix,
};
use ddg_repulsive_surfaces::metric::{FaceOperators, FractionalForm, SobolevWeights};
use ddg_repulsive_surfaces::remesh::{CollapseChecks, EditMesh, FlipChecks, RefusalKind};
use ddg_repulsive_surfaces::solve::{
    DensePreconditioner, DenseSaddle, LinearOperator, SchurComplement,
};
use ddg_repulsive_surfaces::constraints::BarycenterMode;
use ddg_repulsive_surfaces::solve::NewtonSteps;
use ddg_repulsive_surfaces::solve::iterative::Backend;
use ddg_repulsive_surfaces::{
    ArmijoRule, EnergyEval, Exponents, FlowOptions, LineSearch, Metric, RsError, Session,
    StartRule, StepStatus, StopRule, Surface, TrialOutcome,
};
use fixture::{Case, Output};

// ----- tolerances (see the table in the module docs) -----

const ENERGY: f64 = 1e-12;
const GRADIENT: f64 = 1e-12;
const DF: f64 = 1e-12;
const MATRIX: f64 = 1e-12;
const PHI: f64 = 1e-12;
const JACOBIAN: f64 = 1e-12;
const SOLVE: f64 = 1e-10;
const POSITION: f64 = 1e-14;

/// One comparison: `max |got − want| / max(max |want|, floor)` against `tol`.
struct Check {
    what: &'static str,
    got: Vec<f64>,
    want: Vec<f64>,
    tol: f64,
    floor: f64,
}

fn check(what: &'static str, got: Vec<f64>, want: Vec<f64>, tol: f64, floor: f64) -> Check {
    Check { what, got, want, tol, floor }
}

/// What the Rust side did with a case.
#[derive(Default)]
struct Run {
    checks: Vec<Check>,
    /// The Rust rank check dropped rows (D10).
    dropped_rows: bool,
    /// A flip succeeded after a local reordering of the two faces (SPEC §10).
    reordered: bool,
    /// Problems found while building `checks` (structural mismatches).
    problems: Vec<String>,
    /// A special reading of the case (module docs).
    special: Option<String>,
}

// ----- statistics -----

#[derive(Default)]
struct Stats {
    cases: usize,
    by_task: BTreeMap<String, usize>,
    ok_cases: usize,
    error_cases: usize,
    skipped: usize,
    special: Vec<String>,
    comparisons: usize,
    /// per quantity: (largest ratio error / tolerance, largest relative error, comparisons)
    worst: BTreeMap<&'static str, (f64, f64, usize)>,
    /// per task: largest relative error over the task's quantities
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
        self.special.extend(o.special);
        self.comparisons += o.comparisons;
        for (k, (r, e, n)) in o.worst {
            let w = self.worst.entry(k).or_insert((0.0, 0.0, 0));
            w.0 = w.0.max(r);
            w.1 = w.1.max(e);
            w.2 += n;
        }
        self.problems.extend(o.problems);
    }

    fn problem(&mut self, name: &str, msg: String) {
        self.problems.push((name.to_string(), msg));
    }

    fn report(&self) -> String {
        let mut s = format!(
            "{} cases ({} answered ok, {} oracle errors, {} skipped), {} comparisons, {} disagreements\n  by task: {:?}\n",
            self.cases,
            self.ok_cases,
            self.error_cases,
            self.skipped,
            self.comparisons,
            self.problems.len(),
            self.by_task
        );
        s += "  quantity: max relative error | worst ratio error/tolerance | comparisons\n";
        for (k, (r, e, n)) in &self.worst {
            s += &format!("  {k:<10} {e:.2e} | {r:.3} | {n}\n");
        }
        if !self.special.is_empty() {
            let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
            for x in &self.special {
                *counts.entry(x.as_str()).or_default() += 1;
            }
            s += &format!("  special readings: {counts:?}\n");
        }
        s
    }
}

// ----- helpers -----

fn flat3(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn hexes(l: &[String]) -> Vec<f64> {
    l.iter().map(|x| fixture::hex(x)).collect()
}

/// The oracle's error kind against the variant of the Rust error.
fn error_matches(kind: &str, e: &RsError) -> bool {
    let v = format!("{e:?}");
    let v = v.split([' ', '{', '(']).next().unwrap_or("");
    match kind {
        "DegenerateFace" => matches!(v, "DegenerateFace" | "RepeatedVertex"),
        "NonOrientable" => v == "VolumeNeedsOrientation",
        k => v == k,
    }
}

fn exponents(c: &Case) -> Result<Exponents, RsError> {
    let (al, be) = c.exponents.unwrap_or((6.0, 12.0));
    Exponents::new(al, be)
}

fn weights(c: &Case) -> SobolevWeights {
    let (wtp, clo, ww) = c.weights.unwrap_or((1.0, 1.0, 0.0));
    SobolevWeights { tp_weight: wtp, low_order: clo != 0.0, willmore_weight: ww }
}

/// The mesh, validated as the oracle validates it (topology, then positions and face areas).
fn mesh(c: &Case) -> Result<Surface, RsError> {
    let s = Surface::new(c.positions.len(), c.faces.clone())?;
    FaceGeometry::new(&s, &c.positions)?;
    Ok(s)
}

// ----- the tasks -----

fn run_energy(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let e = exponents(c)?;
    let pool = make_pool(Some(1))?;
    let geom = FaceGeometry::new(&s, &c.positions)?;
    let (en, g) = energy_and_differential(&pool, &s, &c.positions, &geom, &e)?;
    let mut r = Run::default();
    r.checks.push(check("energy", vec![en], vec![o.scalar("energy")], ENERGY, 1e-300));
    r.checks.push(check("gradient", flat3(&g), o.vectors("d"), GRADIENT, 1e-300));
    Ok(r)
}

fn run_df(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let ops = FaceOperators::new(&s, &c.positions)?;
    let mut got = Vec::new();
    let mut want = Vec::new();
    for l in o.get("df") {
        let k: usize = l[0].parse().unwrap();
        let f: usize = l[1].parse().unwrap();
        let u = &c.u[k];
        let ub = ops.average(u)[f];
        let g = ops.gradient(u)[f];
        got.extend([ub, g[0], g[1], g[2]]);
        want.extend(hexes(&l[2..6]));
    }
    let mut r = Run::default();
    r.checks.push(check("df", got, want, DF, 1e-300));
    Ok(r)
}

fn run_metric(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let e = exponents(c)?;
    let w = weights(c);
    let pool = make_pool(Some(1))?;
    let dm = dense_metric(&pool, &s, &c.positions, &e, &w, None)?;
    let flat_o = |k: &str| o.matrix(k).into_iter().flatten().collect::<Vec<f64>>();
    let mut r = Run::default();
    r.checks.push(check("A", dm.a.as_slice().to_vec(), flat_o("A"), MATRIX, 1e-300));
    if !c.matrices_a_only {
        let ops = FaceOperators::new(&s, &c.positions)?;
        let tau = c.tau.unwrap_or(0.75);
        let lt = fractional_operator(&pool, &ops, tau)?;
        let wmw = willmore_matrix(&dm.w, &dm.mass)?;
        r.checks.push(check("L", lt.as_slice().to_vec(), flat_o("L"), MATRIX, 1e-300));
        r.checks.push(check("B", dm.b.as_slice().to_vec(), flat_o("B"), MATRIX, 1e-300));
        r.checks.push(check("B0", dm.b0.as_slice().to_vec(), flat_o("B0"), MATRIX, 1e-300));
        r.checks.push(check("W", to_dense(&dm.w).as_slice().to_vec(), flat_o("W"), MATRIX, 1e-300));
        r.checks.push(check(
            "WMW",
            to_dense(&wmw).as_slice().to_vec(),
            flat_o("WMW"),
            MATRIX,
            1e-300,
        ));
        r.checks.push(check("M", dm.mass.clone(), hexes(&o.get("M")[0]), MATRIX, 1e-300));
        r.checks.push(check("s", vec![e.s()], vec![o.scalar("s")], 1e-14, 1e-300));
    }
    Ok(r)
}

/// The constraint set; if only the translation is unfixed (a saddle concern in the oracle), pins on
/// every component are added after the oracle's rows.
fn constraint_set(
    s: &Surface,
    c: &Case,
    pin_if_free: bool,
) -> Result<(ConstraintSet, bool), RsError> {
    match ConstraintSet::new(s, &c.positions, &c.constraints, Default::default()) {
        Ok(set) => Ok((set, false)),
        Err(RsError::TranslationNotFixed { .. }) if pin_if_free => {
            let mut cons = c.constraints.clone();
            for comp in s.components() {
                cons.push(Constraint::Pins { vertices: vec![comp.vertices()[0]], targets: None });
            }
            Ok((ConstraintSet::new(s, &c.positions, &cons, Default::default())?, true))
        }
        Err(e) => Err(e),
    }
}

fn run_constraints(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let (set, _) = constraint_set(&s, c, true)?;
    let mut r = Run { dropped_rows: !set.dropped_rows().is_empty(), ..Run::default() };
    if o.error_kind().is_some() {
        return Ok(r);
    }
    let ncomp: usize = o.get("mesh")[0][2].parse().unwrap();
    if s.n_components() != ncomp {
        r.problems.push(format!("components: oracle {ncomp}, Rust {}", s.n_components()));
    }
    let vcomp: Vec<u32> = o.get("vcomp")[0].iter().map(|x| x.parse().unwrap()).collect();
    if s.vertex_components() != vcomp.as_slice() {
        r.problems.push("component numbering differs".into());
    }
    let rows_out = o.get("row");
    let nr = rows_out.len();
    if set.dropped_rows().iter().any(|&x| x < nr) {
        r.problems.push(format!("Rust dropped rows {:?}", set.dropped_rows()));
    }
    let phi = set.values(&c.positions)?;
    let rows: Vec<usize> = (0..nr).collect();
    let j = set.jacobian(&c.positions, &rows)?;
    let n3 = 3 * s.n_vertices();
    let mut jd = vec![0.0; nr * n3];
    for (row, col, x) in j.triplets() {
        jd[row as usize * n3 + col as usize] += x;
    }
    let mut want_phi = Vec::new();
    for (i, l) in rows_out.iter().enumerate() {
        if set.rows()[i].kind.name() != l[1] {
            r.problems.push(format!("row {i}: kind {} vs oracle {}", set.rows()[i].kind.name(), l[1]));
        }
        if set.rows()[i].kind.is_positional() != (l[2] == "positional") {
            r.problems.push(format!("row {i}: positional/nonlinear differs"));
        }
        want_phi.push(fixture::hex(&l[3]));
    }
    // Φ is a difference (value − target): relative to the scale of the values.
    let scale = set.rows()[..nr].iter().fold(1.0f64, |m, row| m.max(row.target.abs()));
    r.checks.push(check("phi", phi[..nr].to_vec(), want_phi, PHI, scale));
    r.checks.push(check("J", jd, o.matrix("J").into_iter().flatten().collect(), JACOBIAN, 1e-300));
    Ok(r)
}

fn run_saddle(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let e = exponents(c)?;
    let w = weights(c);
    let pool = make_pool(Some(1))?;
    let (set, _) = constraint_set(&s, c, false)?;
    let mut r = Run { dropped_rows: !set.dropped_rows().is_empty(), ..Run::default() };
    if o.error_kind().is_some() {
        return Ok(r);
    }
    let dm = dense_metric(&pool, &s, &c.positions, &e, &w, None)?;
    let geom = FaceGeometry::new(&s, &c.positions)?;
    let (_, g) = energy_and_differential(&pool, &s, &c.positions, &geom, &e)?;
    let g = flat3(&g);
    let q = set.positional_rows();
    let nl = set.nonlinear_rows();
    let rows = &o.get("rows")[0];
    if q.len() != rows[0].parse::<usize>().unwrap() || nl.len() != rows[1].parse::<usize>().unwrap()
    {
        r.problems.push(format!("row split ({}, {}) vs oracle {rows:?}", q.len(), nl.len()));
        return Ok(r);
    }
    let cq = set.jacobian(&c.positions, &q)?;
    let cn = set.jacobian(&c.positions, &nl)?;
    let mut saddle = DenseSaddle::new(&pool, &dm.a, &cq)?;
    let schur = SchurComplement::new(&mut saddle, nl.clone(), cn)?;
    let dir = schur.direction(&mut saddle, &g)?;
    let phi = set.values(&c.positions)?;
    let phin: Vec<f64> = nl.iter().map(|&i| phi[i]).collect();
    let h = schur.correction(&phin);
    let gx: f64 = g.iter().zip(&dir.x).map(|(a, b)| a * b).sum();
    r.checks.push(check("g", g.clone(), o.vectors("g"), GRADIENT, 1e-300));
    r.checks.push(check("x", dir.x.clone(), o.vectors("x"), SOLVE, 1e-300));
    // ⟨g, x⟩ can vanish at critical points: compare relative to ‖g‖‖x‖.
    let gn = g.iter().map(|v| v * v).sum::<f64>().sqrt();
    let xn = dir.x.iter().map(|v| v * v).sum::<f64>().sqrt();
    r.checks.push(check("gx", vec![gx], vec![o.scalar("gx")], SOLVE, gn * xn));
    if !nl.is_empty() {
        r.checks.push(check("phin", phin, hexes(&o.get("phin")[0]), SOLVE, 1.0));
        r.checks.push(check("S", schur.matrix().to_vec(), o.matrix("S").concat(), SOLVE, 1e-300));
        r.checks.push(check("h", h, o.vectors("h"), SOLVE, 1e-300));
    }
    Ok(r)
}

fn run_precond(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let e = exponents(c)?;
    let pool = make_pool(Some(1))?;
    let (set, _) = constraint_set(&s, c, false)?;
    let q = set.positional_rows();
    let mut r = Run { dropped_rows: !set.dropped_rows().is_empty(), ..Run::default() };
    if o.error_kind().is_some() {
        return Ok(r);
    }
    if q.len() != o.get("rows")[0][0].parse::<usize>().unwrap() {
        r.problems.push("positional row count differs".into());
        return Ok(r);
    }
    let cq = set.jacobian(&c.positions, &q)?;
    let prec = DensePreconditioner::new(
        &pool,
        &s,
        &c.positions,
        &e,
        &cq,
        c.eps_rel.unwrap_or(1e-10),
        FractionalForm::FacePairs,
    )?;
    let mut rv = flat3(&c.r);
    rv.extend(vec![0.0; q.len()]);
    let z = prec.apply_full(&rv)?;
    let nd = 3 * s.n_vertices();
    r.checks.push(check("eps", vec![prec.epsilon()], vec![o.scalar("eps")], MATRIX, 1e-300));
    r.checks.push(check(
        "Lfr",
        prec.fractional().as_slice().to_vec(),
        o.matrix("Lfr").concat(),
        MATRIX,
        1e-300,
    ));
    r.checks.push(check("z", z[..nd].to_vec(), o.vectors("z"), SOLVE, 1e-300));
    // L′ itself: the operator of the factored saddle against the printed matrix.
    let lp = o.matrix("Lp");
    let n = lp.len();
    let op = prec.laplacian_saddle().operator();
    let mut got = Vec::with_capacity(n * n);
    let mut want = Vec::with_capacity(n * n);
    let mut ej = vec![0.0; n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        ej.iter_mut().for_each(|x| *x = 0.0);
        ej[j] = 1.0;
        op.apply(&ej, &mut col)?;
        for i in 0..n {
            got.push(col[i]);
            want.push(lp[i][j]);
        }
    }
    r.checks.push(check("Lp", got, want, SOLVE, 1e-300));
    Ok(r)
}

fn remesh_error_matches(kind: &str, k: RefusalKind) -> bool {
    match kind {
        "NonManifoldResult" => matches!(k, RefusalKind::RepeatedFace | RefusalKind::NonManifoldResult),
        "FlipNonManifold" => k == RefusalKind::ExistingEdge,
        other => k.name() == other,
    }
}

/// Remeshing: `Err(Ok(..))` never; the result is a `Run` or the text of a mismatch. A refusal of
/// the oracle is checked here (it is not an `RsError`).
fn run_remesh(c: &Case, o: &Output) -> Result<Run, String> {
    let s = Surface::new(c.positions.len(), c.faces.clone()).map_err(|e| e.to_string())?;
    let (op, a, b) = c.op.clone().expect("op line");
    let mut m = EditMesh::new(&s, &c.positions, &c.pinned).map_err(|e| e.to_string())?;
    let res = match op.as_str() {
        "split" => m.split(a, b).map(|_| false),
        "collapse" => m.collapse(a, b, CollapseChecks { foldover: c.foldover }, None).map(|_| false),
        "flip" => m
            .flip(a, b, FlipChecks { delaunay: c.delaunay, foldover: c.foldover })
            .map(|r| r.reordered),
        other => panic!("unknown op {other}"),
    };
    let mut r = Run::default();
    match (o.error_kind(), res) {
        (Some(kind), Err(refusal)) => {
            if !remesh_error_matches(kind, refusal.kind) {
                r.problems.push(format!("oracle {kind}, Rust {refusal}"));
            }
        }
        (Some(kind), Ok(reordered)) => {
            // SPEC §10: the Rust flip reorders a same-direction pair locally; the oracle refuses.
            if kind != "InconsistentOrientation" || !reordered {
                r.problems.push(format!("oracle {kind}, Rust succeeded (reordered: {reordered})"));
            } else {
                r.reordered = true;
                if let Err(e) = m.finish() {
                    r.problems.push(format!("not valid after a reordered flip: {e}"));
                }
            }
        }
        (None, Err(refusal)) => r.problems.push(format!("oracle ok, Rust refused: {refusal}")),
        (None, Ok(reordered)) => {
            if reordered {
                r.problems.push("Rust reordered, oracle did not need to".into());
                return Ok(r);
            }
            let res = m.finish().map_err(|e| e.to_string())?;
            if res.surface.faces() != o.mesh_faces().as_slice() {
                r.problems.push("faces differ".into());
            }
            let want = o.mesh_positions();
            if res.positions.len() != want.len() {
                r.problems.push("vertex count differs".into());
            } else {
                r.checks.push(check("pos", flat3(&res.positions), flat3(&want), POSITION, 1.0));
            }
            let parents: Vec<u32> = o.get("parent").iter().map(|l| l[1].parse().unwrap()).collect();
            if res.face_parents != parents {
                r.problems.push("face parents differ".into());
            }
            let vmap: Vec<Vec<u32>> = (0..res.vertex_map.n_new() as u32)
                .map(|v| res.vertex_map.entry(v).0.to_vec())
                .collect();
            let want_vmap: Vec<Vec<u32>> = o
                .get("vmap")
                .iter()
                .map(|l| l[1..].iter().map(|x| x.parse().unwrap()).collect())
                .collect();
            if vmap != want_vmap {
                r.problems.push("vertex map differs".into());
            }
        }
    }
    Ok(r)
}

/// One dense projected-gradient step through [`Session`] (module docs for the mapping).
fn run_step(c: &Case, o: &Output) -> Result<Run, RsError> {
    let s = mesh(c)?;
    let e = exponents(c)?;
    let w = weights(c);
    let wtp = w.tp_weight;
    let armijo = c.line("armijo").map_or(1e-2, |l| fixture::hex(&l[1]));
    let start = match c.line("start") {
        Some(l) if l[1] == "abs" => StartRule::Absolute { tau: fixture::hex(&l[2]) / wtp },
        Some(l) if l[1] == "maxdisp" => StartRule::MaxDisplacement { lengths: fixture::hex(&l[2]) },
        Some(l) if l[1] == "reference" => StartRule::Reference,
        _ => StartRule::MaxDisplacement { lengths: 0.1 },
    };
    let max_halvings = c.line("maxhalv").map_or(60, |l| l[1].parse().unwrap());
    let newton = match c.line("newton") {
        Some(l) if l[1] == "fixed" => NewtonSteps::Fixed(l[2].parse().unwrap()),
        Some(l) => NewtonSteps::UntilTolerance { max: l[2].parse().unwrap() },
        None => NewtonSteps::UntilTolerance { max: 5 },
    };
    let options = FlowOptions {
        exponents: e,
        energy_eval: EnergyEval::Exact,
        tp_weight: wtp,
        metric: Metric::Sobolev { low_order: w.low_order },
        backend: Backend::Dense,
        constraints: c.constraints.clone(),
        barycenter_mode: Some(BarycenterMode::PerComponent),
        line_search: LineSearch {
            armijo: ArmijoRule::Standard { c: armijo },
            start,
            tau_min_rel: 1e-12,
            max_halvings: Some(max_halvings),
            ccd: None,
            guard_projection: false,
        },
        newton_steps: newton,
        remesh: None,
        stop: StopRule::Off,
        max_iterations: 1,
        threads: Some(1),
        ..FlowOptions::default()
    };
    let mut session = Session::new(s, c.positions.clone(), options)?;
    let rep = session.step()?;
    let mut r = Run::default();
    if o.error_kind().is_some() {
        return Ok(r);
    }
    let lbar = session.mean_edge_length0();
    r.checks.push(check("F", vec![rep.energy_before / wtp], vec![o.scalar("energy0")], SOLVE, 1e-300));
    r.checks.push(check("slope", vec![rep.sobolev_norm_sq / (wtp * wtp)], vec![o.scalar("slope")], SOLVE, 1e-300));
    r.checks.push(check("xnorm", vec![rep.direction_norm / wtp], vec![o.scalar("xnorm")], SOLVE, 1e-300));
    r.checks.push(check("lbar", vec![lbar], vec![o.scalar("lbar")], MATRIX, 1e-300));
    let trials = o.get("trial");
    let result = &o.get("result")[0];
    let accepted = result[0] == "accepted";
    let halvings: u32 = result[2].parse().unwrap();
    if !trials.is_empty() {
        r.checks.push(check("tau", vec![rep.tau0 * wtp], vec![o.scalar("tau0")], SOLVE, 1e-300));
    }
    // Trials: outcome exactly (borderline Armijo margins excepted), step sizes and energies.
    let n = trials.len().min(rep.trials.len());
    for (i, t) in trials.iter().enumerate().take(n) {
        let got = &rep.trials[i];
        let want_outcome = match t[2].as_str() {
            "accepted" => TrialOutcome::Accepted,
            "armijo-failed" => TrialOutcome::ArmijoFailed,
            _ => TrialOutcome::NonFinite,
        };
        if got.outcome != want_outcome {
            let margin = fixture::hex(&t[4]);
            if margin.abs() < 1e-9 {
                r.special = Some(format!("borderline Armijo (margin {margin:.1e})"));
                return Ok(r);
            }
            r.problems.push(format!("trial {i}: {:?} vs oracle {}", got.outcome, t[2]));
            return Ok(r);
        }
        r.checks.push(check("tau", vec![got.tau * wtp], vec![fixture::hex(&t[1])], SOLVE, 1e-300));
        if want_outcome != TrialOutcome::NonFinite {
            r.checks.push(check("F", vec![got.energy / wtp], vec![fixture::hex(&t[3])], SOLVE, 1e-300));
        }
    }
    if trials.len() != rep.trials.len() {
        r.problems.push(format!("{} trials vs oracle {}", rep.trials.len(), trials.len()));
    }
    if (rep.status == StepStatus::Accepted) != accepted || rep.halvings != halvings {
        r.problems.push(format!(
            "result {:?} after {} halvings vs oracle {} after {halvings}",
            rep.status, rep.halvings, result[0]
        ));
    }
    r.checks.push(check("F", vec![rep.energy_trial / wtp], vec![o.scalar("energy1")], SOLVE, 1e-300));
    let newton: u32 = o.get("newton")[0][0].parse().unwrap();
    if rep.newton_steps != newton {
        r.problems.push(format!("{} Newton rounds vs oracle {newton}", rep.newton_steps));
    }
    let set = ConstraintSet::new(
        session.surface(),
        session.positions(),
        &session.constraints(),
        BarycenterMode::PerComponent,
    )?;
    let phi = set.values(session.positions())?;
    let want_phi = hexes(&o.get("phi")[0]);
    let scale = set.rows().iter().fold(1.0f64, |m, row| m.max(row.target.abs()));
    r.checks.push(check("phi", phi, want_phi, SOLVE, scale));
    r.checks.push(check("F", vec![rep.energy_after / wtp], vec![o.scalar("energy-final")], SOLVE, 1e-300));
    let want_x = o.vectors("x");
    r.checks.push(check("pos", flat3(session.positions()), want_x, SOLVE, 1.0));
    Ok(r)
}

// ----- the replay of one case -----

fn replay_case(name: &str, input: &str, output: &str, stats: &mut Stats) {
    let c = fixture::parse_case(input);
    let o = fixture::parse_output(output);
    stats.cases += 1;
    *stats.by_task.entry(c.task.clone()).or_default() += 1;
    if c.task == "step" && c.weights.is_some_and(|w| w.2 != 0.0) {
        // A Willmore term in the metric without a Willmore energy: not expressible (module docs).
        stats.skipped += 1;
        return;
    }
    let oracle_error = o.error_kind().map(str::to_string);
    if oracle_error.as_deref() == Some("Unsupported") || oracle_error.as_deref() == Some("ValueCount")
    {
        stats.skipped += 1;
        return;
    }
    if o.status != "ok" && oracle_error.is_none() {
        stats.problem(name, format!("unexpected oracle status {:?}", o.status));
        return;
    }
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<Run, String> {
        match c.task.as_str() {
            "remesh-ops" => run_remesh(&c, &o).map_err(|e| format!("Rust setup failed: {e}")),
            task => {
                let r = match task {
                    "energy" => run_energy(&c, &o),
                    "df" => run_df(&c, &o),
                    "metric" => run_metric(&c, &o),
                    "constraints" => run_constraints(&c, &o),
                    "saddle" => run_saddle(&c, &o),
                    "precond" => run_precond(&c, &o),
                    "step" => run_step(&c, &o),
                    other => return Err(format!("unknown task {other}")),
                };
                match (&oracle_error, r) {
                    (Some(kind), Err(e)) => {
                        if error_matches(kind, &e) {
                            Ok(Run::default())
                        } else {
                            Err(format!("oracle error {kind}, Rust error {e}"))
                        }
                    }
                    (Some(kind), Ok(run)) => {
                        if kind == "SingularSaddle" && run.dropped_rows {
                            // D10: the rank check repairs the redundancy the oracle refuses.
                            Ok(Run { dropped_rows: true, ..Run::default() })
                        } else {
                            Err(format!("oracle error {kind}, Rust succeeded"))
                        }
                    }
                    (None, Err(e)) => Err(format!("oracle ok, Rust error: {e}")),
                    (None, Ok(run)) => Ok(run),
                }
            }
        }
    }));
    let run = match result {
        Ok(Ok(run)) => run,
        Ok(Err(msg)) => {
            stats.problem(name, msg);
            return;
        }
        Err(_) => {
            stats.problem(name, "Rust panicked".into());
            return;
        }
    };
    if oracle_error.is_some() {
        stats.error_cases += 1;
    } else {
        stats.ok_cases += 1;
    }
    if run.dropped_rows && oracle_error.is_some() {
        stats.special.push("singular-saddle repaired by the rank check".into());
    }
    if run.reordered {
        stats.special.push("flip reordered by Rust".into());
    }
    if let Some(sp) = &run.special {
        stats.special.push(sp.clone());
    }
    for p in &run.problems {
        stats.problem(name, p.clone());
    }
    for ck in run.checks {
        if ck.got.len() != ck.want.len() {
            stats.problem(name, format!("{}: length {} vs oracle {}", ck.what, ck.got.len(), ck.want.len()));
            continue;
        }
        let scale = ck.want.iter().fold(ck.floor, |m, x| m.max(x.abs()));
        let err = ck.got.iter().zip(&ck.want).fold(0.0f64, |m, (a, b)| {
            let d = (a - b).abs();
            if d.is_nan() || m.is_nan() { f64::NAN } else { m.max(d) }
        }) / scale;
        let ratio = err / ck.tol;
        stats.comparisons += 1;
        let w = stats.worst.entry(ck.what).or_insert((0.0, 0.0, 0));
        w.2 += 1;
        if ratio.is_nan() {
            w.0 = f64::INFINITY;
            w.1 = f64::INFINITY;
        } else {
            w.0 = w.0.max(ratio);
            w.1 = w.1.max(err);
        }
        if ratio.is_nan() || ratio > 1.0 {
            stats.problem(name, format!("{}: relative error {err:.3e} > tolerance {:.0e}", ck.what, ck.tol));
        }
    }
}

fn replay_stream(inputs: &str, outputs: &str) -> Stats {
    let mut stats = Stats::default();
    let outs: BTreeMap<String, String> = fixture::blocks(outputs).into_iter().collect();
    for (name, input) in fixture::blocks(inputs) {
        match outs.get(&name) {
            Some(out) => replay_case(&name, &input, out, &mut stats),
            None => stats.problem(&name, "the oracle gave no answer".into()),
        }
    }
    stats
}

// ----- tests -----

#[test]
fn fixtures_are_reproduced() {
    let mut stats = Stats::default();
    let mut per_task: BTreeMap<String, usize> = BTreeMap::new();
    for (name, task, _) in fixture::manifest() {
        let input = std::fs::read_to_string(fixture::dir().join(format!("{name}.in"))).unwrap();
        let output = std::fs::read_to_string(fixture::dir().join(format!("{name}.out"))).unwrap();
        // Fixture files are single cases without the `end` line semantics of a stream.
        replay_case(&name, &input, &output, &mut stats);
        *per_task.entry(task).or_default() += 1;
    }
    println!("{}", stats.report());
    for (n, m) in &stats.problems {
        println!("DISAGREEMENT {n}: {m}");
    }
    assert!(stats.problems.is_empty(), "{} disagreements", stats.problems.len());
    // Minimum counts (the fixtures of the oracle: 70 cases, 7 of them `step`).
    let replayed = stats.ok_cases + stats.error_cases;
    assert!(replayed >= 67, "only {replayed} fixtures replayed");
    for (task, min) in [
        ("energy", 14),
        ("df", 4),
        ("metric", 7),
        ("constraints", 8),
        ("saddle", 8),
        ("precond", 3),
        ("remesh-ops", 16),
        ("step", 7),
    ] {
        assert!(per_task.get(task).copied().unwrap_or(0) >= min, "too few {task} fixtures");
    }
    assert!(stats.special.iter().any(|s| s.starts_with("singular")), "rs-error-singular");
    assert!(stats.special.iter().any(|s| s.starts_with("flip")), "rs-mobius-flip");
}

/// The `step` fixtures (7) replayed alone: every one answered and reproduced.
#[test]
fn step_fixtures_are_reproduced() {
    let mut stats = Stats::default();
    let mut n = 0;
    for (name, task, _) in fixture::manifest() {
        if task != "step" {
            continue;
        }
        n += 1;
        let input = std::fs::read_to_string(fixture::dir().join(format!("{name}.in"))).unwrap();
        let output = std::fs::read_to_string(fixture::dir().join(format!("{name}.out"))).unwrap();
        replay_case(&name, &input, &output, &mut stats);
    }
    println!("{}", stats.report());
    for (n, m) in &stats.problems {
        println!("DISAGREEMENT {n}: {m}");
    }
    assert_eq!(n, 7);
    assert!(stats.problems.is_empty(), "{} disagreements", stats.problems.len());
    assert_eq!(stats.ok_cases, 7, "every step fixture is answered ok and reproduced");
    assert!(stats.special.is_empty(), "no special readings expected: {:?}", stats.special);
}

fn oracle_path() -> PathBuf {
    match std::env::var("DDG_ORACLE") {
        Ok(p) => PathBuf::from(p),
        Err(_) => PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../lean/oracle/.lake/build/bin/ddg-oracle"),
    }
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn fuzz_chunk(oracle: &Path, seed: usize, first: usize, count: usize) -> (Stats, String) {
    let generated = Command::new(oracle)
        .args(["rs-gen", &seed.to_string(), &count.to_string(), &first.to_string()])
        .output()
        .expect("run rs-gen");
    assert!(generated.status.success(), "rs-gen failed");
    let inputs = String::from_utf8(generated.stdout).unwrap();
    let mut child = Command::new(oracle)
        .arg("rs")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("run rs");
    let mut stdin = child.stdin.take().unwrap();
    let text = inputs.clone();
    let writer = std::thread::spawn(move || stdin.write_all(text.as_bytes()).unwrap());
    let answer = child.wait_with_output().expect("wait for rs");
    writer.join().unwrap();
    assert!(answer.status.success(), "rs failed");
    let outputs = String::from_utf8(answer.stdout).unwrap();
    (replay_stream(&inputs, &outputs), inputs)
}

#[test]
#[ignore = "needs the built oracle (cd lean/oracle && lake build); see the module docs"]
fn fuzz_against_the_oracle() {
    let oracle = oracle_path();
    assert!(oracle.exists(), "{} is not built", oracle.display());
    let cases = env_usize("DDG_RS_FUZZ_CASES", 2000);
    let seed = env_usize("DDG_RS_FUZZ_SEED", 1);
    let threads = env_usize(
        "DDG_RS_FUZZ_THREADS",
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
    );
    let chunk = 10;
    let n_chunks = cases.div_ceil(chunk);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let total = std::sync::Mutex::new(Stats::default());
    let bad_inputs = std::sync::Mutex::new(BTreeMap::<String, String>::new());
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
                        let blocks = fixture::blocks(&inputs);
                        let mut map = bad_inputs.lock().unwrap();
                        for (name, _) in &stats.problems {
                            if let Some((_, text)) = blocks.iter().find(|(n, _)| n == name) {
                                map.insert(name.clone(), format!("{text}end\n"));
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
    for (name, msg) in total.problems.iter().take(40) {
        println!("DISAGREEMENT {name}: {msg}");
    }
    let map = bad_inputs.into_inner().unwrap();
    if !map.is_empty() {
        let dir = std::env::var("DDG_RS_FUZZ_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from(env!("CARGO_TARGET_TMPDIR")));
        std::fs::create_dir_all(&dir).unwrap();
        for (name, text) in &map {
            std::fs::write(dir.join(format!("{name}.in")), text).unwrap();
        }
        println!("inputs of the disagreeing cases written to {}", dir.display());
    }
    assert!(total.problems.is_empty(), "{} disagreements", total.problems.len());
}
