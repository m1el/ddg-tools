//! Differential test against the Lean oracle (`lean/oracle`, command `dev`), SPEC §9 and §10.1.
//!
//! The oracle is an independent Float implementation of the six energies and their gradients, of
//! the first L-BFGS step (scaled steepest descent plus the Lewis-Overton weak Wolfe line search)
//! and of 4-1 subdivision. It was written from `SPEC.md` without this crate's code; this test is
//! written from the oracle's README (section "Developability") and this crate's public API.
//!
//! * `fixtures_are_reproduced` (always runs): replays the committed golden fixtures
//!   `lean/oracle/fixtures/dev/*.in|out`.
//! * `fuzz_against_the_oracle` (`#[ignore]`): asks the built oracle for fresh seeded cases
//!   (`ddg-oracle dev-gen SEED COUNT FIRST`), has the oracle answer them (`ddg-oracle dev`) and
//!   replays every answer, in parallel over all cores. Run
//!
//!   ```text
//!   (cd lean/oracle && lake build)
//!   DDG_DEV_FUZZ_CASES=20000 DDG_DEV_FUZZ_SEED=1 \
//!     cargo test --release -p ddg-developability --test oracle -- --ignored --nocapture
//!   ```
//!
//!   Environment: `DDG_DEV_FUZZ_CASES` (default 4000), `DDG_DEV_FUZZ_SEED` (default 1),
//!   `DDG_DEV_FUZZ_THREADS` (default: all cores), `DDG_ORACLE` (path of the executable, default
//!   `lean/oracle/.lake/build/bin/ddg-oracle`), `DDG_DEV_FUZZ_DIR` (where disagreeing inputs are
//!   written, default the target tmp dir).
//!
//! # Tolerances
//!
//! Everything is compared to rounding, not to a modelling error: both sides evaluate the same
//! formulas of the SPEC in `f64`, with different summation orders and a different eigensolver
//! implementation.
//!
//! * energies: `|a - b| <= 1e-9 max(|a|,|b|) + 1e-13`. The relative part is 7 orders above the
//!   rounding of a well conditioned star (`1e-16`); the absolute part covers near-flat stars,
//!   whose energy is a square of normal deviations of size `eps` and is therefore only known to
//!   about `1e-16 sqrt(E)` plus the `1e-17` absolute error of a 3x3 eigenvalue.
//! * gradients: per vertex `|g - g'| <= 1e-9 G + 1e-12 / mean_edge_length`, `G` the largest vertex
//!   gradient of the mesh. Normwise relative to the largest entry: the gradient of a near-flat
//!   star is `O(sqrt(E))` with rounding errors of the normals `1e-16 / h` that are not relative to
//!   the entry itself. The intrinsic energies amplify rounding by `1/sin(phi)` of the angle between
//!   a face normal and the vertex normal (allowed down to `pi - phi = 1e-6`, SPEC D12), so their
//!   worst observed disagreement is the largest (`~1e-9` of `G` over millions of vertices, still
//!   below the tolerance).
//! * ties (SPEC 3.8): a gradient entry is not compared if the star of the vertex or of one of its
//!   neighbours (those are the stars it depends on) has a margin `<= 1e-9 max(1, e)` on either
//!   side, because there the active piece is chosen by rounding. Energies are compared everywhere.
//! * step: the oracle models SPEC 4.1-4.3 (scaled steepest descent plus the weak Wolfe search), so
//!   the optimizer's admissibility floors (`min_face_angle`, `min_fold_angle`) are switched off;
//!   positions to `1e-9` mean edge lengths, energy as above, the trial count and the kind
//!   of ending exactly; not compared if the start has a near tie anywhere.
//! * subdivision: faces and midpoint ids exactly, positions to `1e-14` of the mesh extent.

// The comparisons are written `!(d <= tol)` on purpose: they must fail for NaN.
#![allow(
    clippy::neg_cmp_op_on_partial_ord,
    clippy::needless_range_loop,
    clippy::collapsible_if
)]

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use ddg_core::TriMesh;
use ddg_developability::{
    DevError, Energy, EnergyOptions, LineSearch, Method, OptimizeOptions, OptimizeStatus, evaluate,
    optimize, subdivide,
};

const ENERGIES: [(Energy, &str); 6] = [
    (Energy::Combinatorial, "combinatorial"),
    (Energy::CombinatorialMax, "combinatorial-max"),
    (Energy::Covariance, "covariance"),
    (Energy::CovarianceMax, "covariance-max"),
    (Energy::Intrinsic, "intrinsic"),
    (Energy::IntrinsicMax, "intrinsic-max"),
];

const ENERGY_REL: f64 = 1e-9;
const ENERGY_ABS: f64 = 1e-13;
const GRADIENT_REL: f64 = 1e-9;
const GRADIENT_ABS: f64 = 1e-12;
const TIE: f64 = 1e-9;
/// Steps from energies or gradients at rounding level are not compared (the direction is the
/// normalized gradient, i.e. amplified noise): `E < NOISE_ENERGY` or `mean_edge_length |g|_inf < NOISE_GRADIENT`.
const NOISE_ENERGY: f64 = 1e-8;
const NOISE_GRADIENT: f64 = 1e-6;
/// Smallest margin of a Wolfe test decision (see the oracle's `trial` lines) that is compared.
const DECISION: f64 = 1e-7;

// ----- the text format -----

type Line = Vec<String>;

fn float(token: &str) -> f64 {
    assert_eq!(token.len(), 16, "a float is 16 hexadecimal digits: {token}");
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|e| panic!("{token}: {e}")))
}

fn hex(x: f64) -> String {
    format!("{:016x}", x.to_bits())
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

struct StepParams {
    initial_step: f64,
    max_displacement: f64,
    c1: f64,
    c2: f64,
    max_doublings: u32,
    max_bisections: u32,
    pinned: Vec<u32>,
}

struct Input {
    name: String,
    task: String,
    energies: Vec<(Energy, &'static str)>,
    include_boundary: bool,
    omit_valence_three: bool,
    step: StepParams,
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
}

fn parse_input(name: &str, lines: &[Line]) -> Input {
    let mut input = Input {
        name: name.to_string(),
        task: String::new(),
        energies: Vec::new(),
        include_boundary: true,
        omit_valence_three: true,
        step: StepParams {
            initial_step: 0.1,
            max_displacement: 0.5,
            c1: 0.0,
            c2: 0.5,
            max_doublings: 30,
            max_bisections: 30,
            pinned: Vec::new(),
        },
        positions: Vec::new(),
        faces: Vec::new(),
    };
    for line in lines {
        match line[0].as_str() {
            "task" => input.task = line[1].clone(),
            "energy" => {
                let e = ENERGIES
                    .iter()
                    .find(|(_, n)| *n == line[1])
                    .unwrap_or_else(|| panic!("unknown energy {}", line[1]));
                input.energies.push(*e);
            }
            "boundary" => input.include_boundary = line[1] == "1",
            "omit3" => input.omit_valence_three = line[1] == "1",
            "step" => {
                input.step.initial_step = float(&line[1]);
                input.step.max_displacement = float(&line[2]);
                input.step.c1 = float(&line[3]);
                input.step.c2 = float(&line[4]);
                input.step.max_doublings = int(&line[5]) as u32;
                input.step.max_bisections = int(&line[6]) as u32;
            }
            "pinned" => input
                .step
                .pinned
                .extend(line[1..].iter().map(|t| int(t) as u32)),
            "p" => input
                .positions
                .push([float(&line[1]), float(&line[2]), float(&line[3])]),
            "f" => input.faces.push([
                int(&line[1]) as u32,
                int(&line[2]) as u32,
                int(&line[3]) as u32,
            ]),
            "mesh" => {}
            other => panic!("{name}: unknown input line {other}"),
        }
    }
    input
}

// ----- statistics -----

#[derive(Default, Clone)]
struct EnergyStats {
    sections: usize,
    errors_agreed: usize,
    vertices: usize,
    /// largest `|a - b|` over all vertices, and over vertices with `|e| > 1e-3` relative
    vertex_abs: f64,
    vertex_rel: f64,
    total_abs: f64,
    total_rel: f64,
    gradient_vertices: usize,
    gradient_skipped: usize,
    /// largest `|g - g'| / G`
    gradient_normwise: f64,
    gradient_abs: f64,
    /// counted vertices with a near tie only here, only in the oracle, on both sides
    ties: [usize; 3],
    worst_gradient_case: String,
}

#[derive(Default, Clone)]
struct Stats {
    cases: usize,
    energy: [EnergyStats; 6],
    steps: usize,
    steps_skipped: usize,
    cap_knife_edge: usize,
    steps_refused: usize,
    steps_by_kind: HashMap<String, usize>,
    step_position: f64,
    step_energy: f64,
    subdivisions: usize,
    subdivision_position: f64,
    invalid_meshes_agreed: usize,
    problems: Vec<(String, String)>,
}

impl Stats {
    fn merge(&mut self, o: Stats) {
        self.cases += o.cases;
        for (a, b) in self.energy.iter_mut().zip(o.energy) {
            a.sections += b.sections;
            a.errors_agreed += b.errors_agreed;
            a.vertices += b.vertices;
            a.vertex_abs = a.vertex_abs.max(b.vertex_abs);
            a.vertex_rel = a.vertex_rel.max(b.vertex_rel);
            a.total_abs = a.total_abs.max(b.total_abs);
            a.total_rel = a.total_rel.max(b.total_rel);
            a.gradient_vertices += b.gradient_vertices;
            a.gradient_skipped += b.gradient_skipped;
            a.gradient_normwise = a.gradient_normwise.max(b.gradient_normwise);
            a.gradient_abs = a.gradient_abs.max(b.gradient_abs);
            for i in 0..3 {
                a.ties[i] += b.ties[i];
            }
            if b.gradient_normwise >= a.gradient_normwise {
                a.worst_gradient_case = b.worst_gradient_case.clone();
            }
        }
        self.steps += o.steps;
        self.steps_skipped += o.steps_skipped;
        self.cap_knife_edge += o.cap_knife_edge;
        self.steps_refused += o.steps_refused;
        for (k, v) in o.steps_by_kind {
            *self.steps_by_kind.entry(k).or_default() += v;
        }
        self.step_position = self.step_position.max(o.step_position);
        self.step_energy = self.step_energy.max(o.step_energy);
        self.subdivisions += o.subdivisions;
        self.subdivision_position = self.subdivision_position.max(o.subdivision_position);
        self.invalid_meshes_agreed += o.invalid_meshes_agreed;
        self.problems.extend(o.problems);
    }

    fn report(&self) -> String {
        let mut s = format!("{} cases replayed\n", self.cases);
        s += "energy                 sections  vertices  max|de| (all)  max rel de (e>1e-3)  max|dE|     grad vertices  skipped(ties)  max|dg|/G   max|dg|\n";
        for (i, (_, name)) in ENERGIES.iter().enumerate() {
            let e = &self.energy[i];
            s += &format!(
                "{name:<22} {:>8} {:>9} {:>14.2e} {:>20.2e} {:>10.2e} {:>14} {:>14} {:>10.2e} {:>10.2e}\n",
                e.sections,
                e.vertices,
                e.vertex_abs,
                e.vertex_rel,
                e.total_abs,
                e.gradient_vertices,
                e.gradient_skipped,
                e.gradient_normwise,
                e.gradient_abs
            );
        }
        s += "near ties among counted vertices (margin <= 1e-9 max(1,e)): only here / only in the oracle / both\n";
        for (i, (_, name)) in ENERGIES.iter().enumerate() {
            let t = self.energy[i].ties;
            s += &format!("  {name:<22} {:>8} {:>8} {:>8}\n", t[0], t[1], t[2]);
        }
        s += "worst gradient disagreements (case, vertex):\n";
        for (i, (_, name)) in ENERGIES.iter().enumerate() {
            s += &format!("  {name:<22} {}\n", self.energy[i].worst_gradient_case);
        }
        s += &format!(
            "steps: {} compared ({} skipped: ties, rounding-level energy or gradient, close decisions; {} capped steps with one trial more or less at the cap), endings {:?}, largest position difference {:.2e} mean edge lengths, energy {:.2e}\n",
            self.steps,
            self.steps_skipped,
            self.cap_knife_edge,
            self.steps_by_kind,
            self.step_position,
            self.step_energy
        );
        s += &format!(
            "subdivisions: {} compared, largest position difference {:.2e}; invalid meshes rejected by both: {}; steps refused by both (undefined energy at the start): {}\n",
            self.subdivisions,
            self.subdivision_position,
            self.invalid_meshes_agreed,
            self.steps_refused
        );
        s += &format!("{} disagreements\n", self.problems.len());
        s
    }
}

// ----- replay of one case -----

fn mean_edge_length(mesh: &TriMesh, positions: &[[f64; 3]]) -> f64 {
    let sum: f64 = mesh
        .edges()
        .iter()
        .map(|e| dist(positions[e[0] as usize], positions[e[1] as usize]))
        .sum();
    sum / mesh.n_edges() as f64
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn replay(input: &Input, out: &[Line], stats: &mut Stats) {
    stats.cases += 1;
    let status = out.iter().find(|l| l[0] == "status").expect("status line");
    let mesh = TriMesh::new(input.positions.len(), input.faces.clone());
    if status[1] != "ok"
        && input.task == "step"
        && (status[2] == "undefined" || status[2] == "degenerate-face")
    {
        // the energy is undefined (or a face degenerate) at the start of the step: `optimize` must refuse
        let Ok(m) = &mesh else {
            stats.problems.push((
                input.name.clone(),
                "TriMesh::new rejects a mesh the oracle accepts".into(),
            ));
            return;
        };
        let (energy, _) = input.energies[0];
        let eopts = EnergyOptions {
            energy,
            include_boundary: input.include_boundary,
            omit_valence_three: input.omit_valence_three,
            threads: None,
        };
        let opts = OptimizeOptions {
            max_iterations: 1,
            gradient_tolerance: 0.0,
            min_face_angle: 0.0,
            min_fold_angle: 0.0,
            ..OptimizeOptions::default()
        };
        let ok = match (
            optimize(m, &input.positions, &eopts, &opts, None),
            status[2].as_str(),
        ) {
            (Err(DevError::UndefinedEnergy { .. }), "undefined") => true,
            (Err(DevError::DegenerateFace { .. }), "degenerate-face") => true,
            _ => false,
        };
        if ok {
            stats.steps_refused += 1;
        } else {
            stats.problems.push((
                input.name.clone(),
                format!(
                    "oracle: step refused ({}), here optimize does not fail the same way",
                    status[2..].join(" ")
                ),
            ));
        }
        return;
    }
    if status[1] != "ok" {
        match mesh {
            Err(_) => stats.invalid_meshes_agreed += 1,
            Ok(_) => stats.problems.push((
                input.name.clone(),
                format!(
                    "the oracle rejects the mesh ({}) but TriMesh::new accepts it",
                    status[2..].join(" ")
                ),
            )),
        }
        return;
    }
    let mesh = match mesh {
        Ok(m) => m,
        Err(e) => {
            stats.problems.push((
                input.name.clone(),
                format!("TriMesh::new rejects the mesh ({e}) but the oracle accepts it"),
            ));
            return;
        }
    };
    match input.task.as_str() {
        "energy" => replay_energy(input, &mesh, out, stats),
        "step" => replay_step(input, &mesh, out, stats),
        "subdivide" => replay_subdivide(input, &mesh, out, stats),
        t => stats
            .problems
            .push((input.name.clone(), format!("unknown task {t}"))),
    }
}

fn replay_energy(input: &Input, mesh: &TriMesh, out: &[Line], stats: &mut Stats) {
    let n = mesh.n_vertices();
    let ell = mean_edge_length(mesh, &input.positions);
    // neighbours (the gradient at p depends on the stars of p and its neighbours)
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for e in mesh.edges() {
        neighbours[e[0] as usize].push(e[1] as usize);
        neighbours[e[1] as usize].push(e[0] as usize);
    }
    let mut at = out
        .iter()
        .position(|l| l[0] == "section")
        .expect("a section");
    for (energy, name) in &input.energies {
        let k = ENERGIES.iter().position(|(_, nm)| nm == name).unwrap();
        let es = &mut stats.energy[k];
        es.sections += 1;
        assert_eq!(out[at][0], "section", "{}: section expected", input.name);
        assert_eq!(out[at][1], *name);
        let result = &out[at + 1];
        assert_eq!(result[0], "result");
        let opts = EnergyOptions {
            energy: *energy,
            include_boundary: input.include_boundary,
            omit_valence_three: input.omit_valence_three,
            threads: None,
        };
        let mine = evaluate(mesh, &input.positions, &opts);
        let problem = |stats: &mut Stats, msg: String| {
            stats
                .problems
                .push((format!("{} [{name}]", input.name), msg))
        };
        if result[1] == "error" {
            at += 2;
            let want = result[2].as_str();
            let agree = match (&mine, want) {
                (Err(DevError::DegenerateFace { .. }), "degenerate-face") => true,
                (Err(DevError::UndefinedEnergy { vertex, .. }), "undefined") => {
                    // which vertex is reported first is not specified; only the kind is compared
                    let _ = (vertex, &result[3]);
                    true
                }
                _ => false,
            };
            if agree {
                stats.energy[k].errors_agreed += 1;
            } else {
                problem(
                    stats,
                    format!(
                        "oracle error '{}', here {:?}",
                        result[2..].join(" "),
                        mine.as_ref().map(|e| e.energy)
                    ),
                );
            }
            continue;
        }
        let ev = match mine {
            Ok(ev) => ev,
            Err(e) => {
                // skip the rest of the section's lines
                at += 2;
                while at < out.len() && out[at][0] != "section" {
                    at += 1;
                }
                problem(stats, format!("oracle ok, here error: {e}"));
                continue;
            }
        };
        let total = float(&result[2]);
        at += 2;
        let es = &mut stats.energy[k];
        let d = (total - ev.energy).abs();
        es.total_abs = es.total_abs.max(d);
        if total.abs().max(ev.energy.abs()) > 1e-3 {
            es.total_rel = es.total_rel.max(d / total.abs().max(ev.energy.abs()));
        }
        let mut local: Vec<String> = Vec::new();
        if !(d <= ENERGY_REL * total.abs().max(ev.energy.abs()) + ENERGY_ABS) {
            local.push(format!(
                "total energy: here {:e}, oracle {:e}",
                ev.energy, total
            ));
        }
        // per-vertex lines: v, g, ch
        let mut o_energy = vec![0.0; n];
        let mut o_margin = vec![0.0; n];
        let mut o_grad = vec![[0.0; 3]; n];
        let mut o_counted = vec![false; n];
        for _ in 0..n {
            let v = &out[at];
            assert_eq!(v[0], "v");
            let i = int(&v[1]);
            o_counted[i] = v[2] == "1";
            o_energy[i] = float(&v[5]);
            o_margin[i] = float(&v[6]);
            let g = &out[at + 1];
            assert_eq!(g[0], "g");
            o_grad[int(&g[1])] = [float(&g[2]), float(&g[3]), float(&g[4])];
            at += 3;
        }
        let es = &mut stats.energy[k];
        let mut gscale = 0.0f64;
        for i in 0..n {
            if o_counted[i] != ev.counted[i] {
                local.push(format!(
                    "vertex {i}: counted here {}, oracle {}",
                    ev.counted[i], o_counted[i]
                ));
            }
            es.vertices += 1;
            let (a, b) = (ev.vertex_energy[i], o_energy[i]);
            if a.is_nan() && b.is_nan() {
                continue;
            }
            let d = (a - b).abs();
            es.vertex_abs = es.vertex_abs.max(d);
            if a.abs().max(b.abs()) > 1e-3 {
                es.vertex_rel = es.vertex_rel.max(d / a.abs().max(b.abs()));
            }
            if !(d <= ENERGY_REL * a.abs().max(b.abs()) + ENERGY_ABS) {
                local.push(format!("vertex {i}: energy here {a:e}, oracle {b:e}"));
            }
            gscale = gscale.max(norm(ev.gradient[i])).max(norm(o_grad[i]));
        }
        for i in 0..n {
            if ev.counted[i] {
                let scale = TIE * 1.0f64.max(ev.vertex_energy[i].abs());
                let (r, o) = (!(ev.vertex_margin[i] > scale), !(o_margin[i] > scale));
                es.ties[if r && o {
                    2
                } else if r {
                    0
                } else if o {
                    1
                } else {
                    3
                }
                .min(2)] += (r || o) as usize;
            }
        }
        // tie-free vertices: the star and all neighbour stars
        let tie = |u: usize| -> bool {
            if !ev.counted[u] {
                return false;
            }
            let scale = 1.0f64.max(ev.vertex_energy[u].abs());
            let m = ev.vertex_margin[u].min(o_margin[u]);
            !(m > TIE * scale) // also true for NaN
        };
        for p in 0..n {
            if tie(p) || neighbours[p].iter().any(|&u| tie(u)) {
                es.gradient_skipped += 1;
                continue;
            }
            es.gradient_vertices += 1;
            let g = ev.gradient[p];
            let h = o_grad[p];
            let d = norm([g[0] - h[0], g[1] - h[1], g[2] - h[2]]);
            es.gradient_abs = es.gradient_abs.max(d);
            if gscale > 1e-6 && d / gscale > es.gradient_normwise {
                es.gradient_normwise = d / gscale;
                es.worst_gradient_case = format!("{} vertex {p}", input.name);
            }
            if !(d <= GRADIENT_REL * gscale + GRADIENT_ABS / ell) {
                local.push(format!(
                    "vertex {p}: gradient here {g:?}, oracle {h:?} (G = {gscale:e})"
                ));
            }
        }
        if !local.is_empty() {
            let more = local.len().saturating_sub(4);
            local.truncate(4);
            if more > 0 {
                local.push(format!("... and {more} more"));
            }
            problem(stats, local.join("; "));
        }
    }
}

fn replay_step(input: &Input, mesh: &TriMesh, out: &[Line], stats: &mut Stats) {
    let (energy, name) = input.energies[0];
    let k = ENERGIES.iter().position(|(_, nm)| *nm == name).unwrap();
    let key = |tag: &str| {
        out.iter()
            .find(|l| l[0] == tag)
            .unwrap_or_else(|| panic!("{}: no {tag}", input.name))
    };
    let result = key("result");
    let o_margin = float(&key("margin")[1]);
    let o_lbar = float(&key("lbar")[1]);
    let eopts = EnergyOptions {
        energy,
        include_boundary: input.include_boundary,
        omit_valence_three: input.omit_valence_three,
        threads: None,
    };
    let _ = k;
    let ev = match evaluate(mesh, &input.positions, &eopts) {
        Ok(ev) => ev,
        Err(_) => {
            stats.problems.push((
                input.name.clone(),
                "oracle ok at the start, here error".into(),
            ));
            return;
        }
    };
    let ell = mean_edge_length(mesh, &input.positions);
    let mut local = Vec::new();
    if (ell - o_lbar).abs() > 1e-12 * ell {
        local.push(format!("mean edge length: here {ell:e}, oracle {o_lbar:e}"));
    }
    let tie = (0..mesh.n_vertices()).any(|v| {
        ev.counted[v] && !(ev.vertex_margin[v] > TIE * 1.0f64.max(ev.vertex_energy[v].abs()))
    }) || !(o_margin > TIE);
    // the decisions of the line search compare energies and slopes with rounding errors of the
    // energy: not meaningful when the energy or the gradient is at rounding level, or when a
    // decision is closer than `DECISION` (relative to the first order change `|g.d| t`)
    let e0 = float(&key("energy0")[1]);
    let gmax = float(&key("gmax")[1]);
    let noise = e0 < NOISE_ENERGY || ell * gmax < NOISE_GRADIENT;
    let close_decision = out.iter().filter(|l| l[0] == "trial").any(|l| {
        let (sm, cm) = (float(&l[6]), float(&l[7]));
        sm.abs() < DECISION || (l[5] != "-1" && cm.abs() < DECISION)
    });
    if tie || noise || close_decision {
        stats.steps_skipped += 1;
        if !local.is_empty() {
            stats.problems.push((input.name.clone(), local.join("; ")));
        }
        return;
    }
    let s = &input.step;
    let opts = OptimizeOptions {
        method: Method::Lbfgs { memory: 10 },
        line_search: LineSearch::WeakWolfe {
            c1: s.c1,
            c2: s.c2,
            max_doublings: s.max_doublings,
            max_bisections: s.max_bisections,
        },
        max_iterations: 1,
        gradient_tolerance: 0.0,
        initial_step: s.initial_step,
        max_displacement: s.max_displacement,
        pinned: s.pinned.clone(),
        // the oracle implements SPEC §4.1-§4.3; the admissibility floors (D31, D32) are an addition of
        // the optimizer that the oracle does not model, so they are switched off here
        min_face_angle: 0.0,
        min_fold_angle: 0.0,
        ..OptimizeOptions::default()
    };
    let res = match optimize(mesh, &input.positions, &eopts, &opts, None) {
        Ok(r) => r,
        Err(e) => {
            stats
                .problems
                .push((input.name.clone(), format!("optimize failed: {e}")));
            return;
        }
    };
    stats.steps += 1;
    let kind = result[1].as_str();
    *stats.steps_by_kind.entry(kind.to_string()).or_default() += 1;
    let o_trials = int(&result[3]);
    let r = &res.report;
    // ending
    match kind {
        "stationary" => {
            if r.iterations != 0 {
                local.push(format!(
                    "oracle: stationary, here {} iterations",
                    r.iterations
                ));
            }
        }
        "failed" => {
            if r.status != OptimizeStatus::LineSearchFailed {
                local.push(format!("oracle: line search failed, here {:?}", r.status));
            }
        }
        _ => {
            let (capped, weak) = (r.capped_steps, r.weak_steps);
            let ok = match kind {
                "wolfe" => capped == 0 && weak == 0,
                "capped" => capped == 1 && weak == 0,
                "weak" => capped == 0 && weak == 1,
                _ => false,
            };
            if !ok || r.iterations != 1 {
                local.push(format!("oracle ending '{kind}', here iterations {} capped {capped} weak {weak} status {:?}", r.iterations, r.status));
            }
            if r.line_search_trials != o_trials {
                // Knife edge: with `initial_step = 0.1` and `max_displacement = 0.2` the cap is
                // `t_max = 2` up to the last bit, and the doubling phase reaches exactly 2; whether
                // that trial is the capped one or one more is evaluated (`alpha >= t_max`) depends
                // on the last bit of `t_max`, i.e. on the rounding of `|d|_inf`. The accepted step
                // is the same to 1e-16, so one trial more or less is accepted for capped steps.
                if kind == "capped" && r.line_search_trials.abs_diff(o_trials) == 1 {
                    stats.cap_knife_edge += 1;
                } else {
                    local.push(format!(
                        "trials: here {}, oracle {o_trials}",
                        r.line_search_trials
                    ));
                }
            }
            let e1 = float(&result[4]);
            let de = (e1 - r.energy_final).abs();
            stats.step_energy = stats.step_energy.max(de);
            if !(de <= ENERGY_REL * e1.abs().max(r.energy_final.abs()) + ENERGY_ABS) {
                local.push(format!(
                    "energy after the step: here {:e}, oracle {e1:e}",
                    r.energy_final
                ));
            }
        }
    }
    let mut worst = 0.0f64;
    for l in out.iter().filter(|l| l[0] == "x") {
        let i = int(&l[1]);
        let p = [float(&l[2]), float(&l[3]), float(&l[4])];
        worst = worst.max(dist(p, res.positions[i]) / ell);
    }
    stats.step_position = stats.step_position.max(worst);
    if !(worst <= 1e-9) {
        local.push(format!(
            "positions after the step differ by {worst:e} mean edge lengths"
        ));
    }
    if !local.is_empty() {
        stats
            .problems
            .push((format!("{} [step {name}]", input.name), local.join("; ")));
    }
}

fn replay_subdivide(input: &Input, mesh: &TriMesh, out: &[Line], stats: &mut Stats) {
    let sub = match subdivide(mesh, &input.positions) {
        Ok(s) => s,
        Err(e) => {
            stats
                .problems
                .push((input.name.clone(), format!("subdivide failed: {e}")));
            return;
        }
    };
    stats.subdivisions += 1;
    let mut local = Vec::new();
    let nv = mesh.n_vertices();
    let header = out.iter().find(|l| l[0] == "mesh").unwrap();
    if int(&header[1]) != sub.positions.len() || int(&header[2]) != sub.mesh.n_faces() {
        local.push(format!(
            "sizes: here {} vertices {} faces, oracle {} {}",
            sub.positions.len(),
            sub.mesh.n_faces(),
            header[1],
            header[2]
        ));
    }
    let extent = input
        .positions
        .iter()
        .flat_map(|p| p.iter().map(|x| x.abs()))
        .fold(0.0f64, f64::max);
    let mut worst = 0.0f64;
    for l in out.iter().filter(|l| l[0] == "x") {
        let i = int(&l[1]);
        if i >= sub.positions.len() {
            local.push(format!("oracle vertex {i} does not exist here"));
            continue;
        }
        let p = [float(&l[2]), float(&l[3]), float(&l[4])];
        worst = worst.max(dist(p, sub.positions[i]));
    }
    stats.subdivision_position = stats.subdivision_position.max(worst / extent.max(1e-300));
    if worst > 1e-14 * extent {
        local.push(format!("positions differ by {worst:e}"));
    }
    let faces: Vec<[u32; 3]> = out
        .iter()
        .filter(|l| l[0] == "f")
        .map(|l| [int(&l[1]) as u32, int(&l[2]) as u32, int(&l[3]) as u32])
        .collect();
    if faces != sub.mesh.faces() {
        local.push("faces differ".into());
    }
    for l in out.iter().filter(|l| l[0] == "em") {
        let (a, b, id) = (int(&l[1]) as u32, int(&l[2]) as u32, int(&l[3]));
        match mesh.find_edge(a as usize, b as usize) {
            Some(e) if sub.edge_midpoint[e] as usize == id => {}
            other => local.push(format!(
                "midpoint of edge ({a}, {b}): oracle {id}, here {:?} (edge {other:?}, nv {nv})",
                other.map(|e| sub.edge_midpoint[e])
            )),
        }
    }
    if !local.is_empty() {
        stats.problems.push((input.name.clone(), local.join("; ")));
    }
}

fn replay_text(inputs: &str, outputs: &str) -> Stats {
    let mut stats = Stats::default();
    let ins = blocks(inputs);
    let outs: HashMap<String, Vec<Line>> = blocks(outputs).into_iter().collect();
    for (name, lines) in ins {
        let input = parse_input(&name, &lines);
        match outs.get(&name) {
            Some(out) => replay(&input, out, &mut stats),
            None => stats
                .problems
                .push((name, "no answer from the oracle".into())),
        }
    }
    stats
}

// ----- fixtures -----

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/dev")
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
    assert!(
        names.len() >= 40,
        "only {} fixtures in {}",
        names.len(),
        dir.display()
    );
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
    assert!(
        total.problems.is_empty(),
        "{} disagreements with the golden fixtures",
        total.problems.len()
    );
}

// ----- fuzzing -----

fn oracle_path() -> PathBuf {
    match std::env::var("DDG_ORACLE") {
        Ok(p) => PathBuf::from(p),
        Err(_) => PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../lean/oracle/.lake/build/bin/ddg-oracle"),
    }
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Generate `count` cases starting at index `first`, answer them, replay them.
fn fuzz_chunk(oracle: &Path, seed: usize, first: usize, count: usize) -> (Stats, String) {
    let generated = Command::new(oracle)
        .args([
            "dev-gen",
            &seed.to_string(),
            &count.to_string(),
            &first.to_string(),
        ])
        .output()
        .expect("run dev-gen");
    assert!(generated.status.success(), "dev-gen failed");
    let inputs = String::from_utf8(generated.stdout).unwrap();
    let mut child = Command::new(oracle)
        .arg("dev")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("run dev");
    let mut stdin = child.stdin.take().unwrap();
    let text = inputs.clone();
    let writer = std::thread::spawn(move || {
        stdin.write_all(text.as_bytes()).unwrap();
    });
    let answer = child.wait_with_output().expect("wait for dev");
    writer.join().unwrap();
    assert!(answer.status.success(), "dev failed");
    let outputs = String::from_utf8(answer.stdout).unwrap();
    (replay_text(&inputs, &outputs), inputs)
}

#[test]
#[ignore = "needs the built oracle (cd lean/oracle && lake build); see the module docs"]
fn fuzz_against_the_oracle() {
    let oracle = oracle_path();
    assert!(oracle.exists(), "{} is not built", oracle.display());
    let cases = env_usize("DDG_DEV_FUZZ_CASES", 4000);
    let seed = env_usize("DDG_DEV_FUZZ_SEED", 1);
    let threads = env_usize(
        "DDG_DEV_FUZZ_THREADS",
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4),
    );
    let chunk = 25;
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
                            let base = name.split(' ').next().unwrap().to_string();
                            if let Some((_, lines)) =
                                blocks(&inputs).into_iter().find(|(n, _)| *n == base)
                            {
                                let mut text = format!("case {base}\n");
                                for l in lines {
                                    text += &(l.join(" ") + "\n");
                                }
                                text += "end\n";
                                map.insert(base, text);
                            }
                        }
                    }
                    total.lock().unwrap().merge(stats);
                }
            });
        }
    });
    let total = total.into_inner().unwrap();
    println!(
        "seed {seed}, {cases} cases, {threads} threads, {:.1} s",
        started.elapsed().as_secs_f64()
    );
    println!("{}", total.report());
    let dir = std::env::var("DDG_DEV_FUZZ_DIR")
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
        println!(
            "inputs of the disagreeing cases written to {}",
            dir.display()
        );
    }
    let _ = hex(0.0);
    assert!(
        total.problems.is_empty(),
        "{} disagreements",
        total.problems.len()
    );
}

// ----- a finding that the comparison above cannot show (both sides share it) -----

/// A boundary vertex with two faces has `e_v = 0` for every energy (SPEC 3.1), so including the
/// boundary in the domain must not change the energy or the gradient of the interior vertices.
/// For `E_lambda` and the intrinsic width it does, when the two faces of a boundary vertex are
/// nearly coplanar: the smallest eigenvalue of its covariance is double up to the square of the
/// normal angle `delta`, the computed eigenvector is an arbitrary vector of the null plane, and
/// the sensitivities `2 theta (x.N) x` are of the size `delta` instead of zero (the oracle has the
/// same noise, which is why it reports the real eigenvalue gap for such stars and the comparison
/// skips them). Run with `--ignored --nocapture`; it prints the relative change of the gradient
/// when the (zero) boundary energies are added, for hinge stars with noise `eps`.
#[test]
#[ignore = "informational: prints the effect; see the comment"]
fn nearly_coplanar_two_face_boundary_stars_perturb_the_gradient() {
    // hinge star: apex 0, ring on the planes z = +- a |x|, two ring points on the crease
    let ring_angles = [90.0f64, 130.0, 180.0, 230.0, 270.0, 320.0, 0.0, 40.0];
    for eps in [1e-12, 1e-9, 1e-6, 1e-3] {
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut noise = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2.0 * eps
        };
        let mut positions = vec![[noise(), noise(), noise()]];
        for a in ring_angles {
            let t = a.to_radians();
            let (x, y) = (1.1 * t.cos(), 1.1 * t.sin());
            let x = if x.abs() < 1e-12 { 0.0 } else { x };
            positions.push([x + noise(), y + noise(), 0.8 * x.abs() + noise()]);
        }
        let k = ring_angles.len();
        let faces: Vec<[u32; 3]> = (0..k)
            .map(|m| [0, 1 + m as u32, 1 + ((m + 1) % k) as u32])
            .collect();
        let mesh = TriMesh::new(positions.len(), faces).unwrap();
        for energy in [Energy::Covariance, Energy::Intrinsic, Energy::Combinatorial] {
            let with = |include_boundary| {
                let opts = EnergyOptions {
                    energy,
                    include_boundary,
                    omit_valence_three: true,
                    threads: None,
                };
                evaluate(&mesh, &positions, &opts).unwrap()
            };
            let (a, b) = (with(true), with(false));
            let g0 = norm(b.gradient[0]).max(1e-300);
            let diff = (0..positions.len())
                .map(|v| {
                    norm([
                        a.gradient[v][0] - b.gradient[v][0],
                        a.gradient[v][1] - b.gradient[v][1],
                        a.gradient[v][2] - b.gradient[v][2],
                    ])
                })
                .fold(0.0f64, f64::max);
            println!(
                "eps {eps:e} {energy:?}: |g| at the apex {g0:.3e}, change when the boundary is included {diff:.3e} (relative {:.3e})",
                diff / g0
            );
        }
    }
}
