//! The developability flow (SPEC §5.1; paper Sec. 4.2 and 5.2): rounds of optimization with
//! optional 4-1 subdivision between them and small-angle remeshing during them.
//!
//! For each [`Round`]: subdivide if asked (pinned ids follow SPEC §5.2), remesh if enabled
//! (SPEC §5.3), then [`optimize`] with the remeshing check after every
//! `check_interval` accepted iterations. When the check finds triangles with a corner angle below
//! the small angle that were not already known to be unfixable, the optimizer is stopped from its
//! progress callback, [`remesh`] is applied, and the optimizer is restarted (its L-BFGS memory is
//! fresh per call) with the remaining iteration budget of the round. A remesh that changes
//! nothing does not restart anything.
//!
//! Nothing is modified in place; the flow returns the final mesh, positions, and for every final
//! vertex the id of the input vertex it descends from.

use std::sync::mpsc;
use std::time::Instant;

use ddg_core::topology::{NONE, TriMesh};

use crate::energy::{EnergyOptions, evaluate_energy};
use crate::error::DevError;
use crate::geometry::{corner_angles, v3};
use crate::optimize::{START_SLACK, check_admissible, 
    Control, OptimizeOptions, OptimizeReport, OptimizeStatus, Progress, optimize,
};
use crate::remesh::{RemeshOptions, RemeshReport, remesh};
use crate::subdivide::subdivide;

/// One round of the flow: an energy, an iteration budget, and whether to subdivide first.
#[derive(Debug, Clone, PartialEq)]
pub struct Round {
    /// The energy minimized in this round (the paper uses `E_P` coarse and `E_λ` fine for
    /// Figs. 23–24).
    pub energy: EnergyOptions,
    /// Iteration budget of the round, over all restarts after remeshing. Overrides
    /// `FlowOptions::optimize.max_iterations`.
    pub max_iterations: usize,
    /// Apply one 4-1 subdivision before the round.
    pub subdivide_before: bool,
}

/// Options of [`flow`] (SPEC §5.1).
#[derive(Debug, Clone, PartialEq)]
pub struct FlowOptions {
    /// The rounds, at least one.
    pub rounds: Vec<Round>,
    /// Optimizer options shared by all rounds. `max_iterations` is replaced by the round's, and
    /// `pinned` must be empty (use [`FlowOptions::pinned`]).
    pub optimize: OptimizeOptions,
    /// Small-angle remeshing; `None` switches it off.
    pub remesh: Option<RemeshOptions>,
    /// Vertices (input ids) that do not move. A midpoint created by subdivision is pinned iff
    /// both ends are; a collapse keeps the pinned vertex.
    pub pinned: Vec<u32>,
    /// `Some(n)`: run the whole flow in a local thread pool of `n` threads (the energies of the
    /// rounds then use it; their own `threads` fields are overridden). `None`: the energies'
    /// own `threads` fields apply. Results do not depend on the thread count (SPEC §3.9).
    pub threads: Option<usize>,
}

impl FlowOptions {
    /// `subdivisions + 1` equal rounds of `iterations` iterations each with `energy`; every round
    /// but the first starts with a 4-1 subdivision. Remeshing on with defaults; nothing pinned.
    pub fn uniform(energy: EnergyOptions, subdivisions: usize, iterations: usize) -> Self {
        FlowOptions {
            rounds: (0..=subdivisions)
                .map(|r| Round { energy: energy.clone(), max_iterations: iterations, subdivide_before: r > 0 })
                .collect(),
            optimize: OptimizeOptions::default(),
            remesh: Some(RemeshOptions::default()),
            pinned: Vec::new(),
            threads: None,
        }
    }

    /// Checks the options.
    ///
    /// # Errors
    /// [`DevError::InvalidOption`] naming the first violation.
    pub fn validate(&self) -> Result<(), DevError> {
        if self.rounds.is_empty() {
            return Err(DevError::InvalidOption { name: "rounds", value: 0.0, requirement: ">= 1 round" });
        }
        for r in &self.rounds {
            r.energy.validate()?;
        }
        self.optimize.validate().map_err(|e| DevError::InvalidOption {
            name: e.name,
            value: e.value,
            requirement: e.requirement,
        })?;
        if !self.optimize.pinned.is_empty() {
            return Err(DevError::InvalidOption {
                name: "optimize.pinned",
                value: self.optimize.pinned.len() as f64,
                requirement: "empty (use FlowOptions::pinned; the flow tracks ids across subdivision and remeshing)",
            });
        }
        if let Some(r) = &self.remesh {
            r.validate()?;
        }
        if self.threads == Some(0) {
            return Err(DevError::InvalidOption { name: "threads", value: 0.0, requirement: ">= 1 (or None)" });
        }
        Ok(())
    }
}

impl Default for FlowOptions {
    /// Two subdivisions (three rounds) of 500 iterations with `E_λ`, remeshing on (≈ the setup
    /// of the paper's Fig. 1).
    fn default() -> Self {
        FlowOptions::uniform(EnergyOptions::default(), 2, 500)
    }
}

/// What happened in one round.
#[derive(Debug, Clone, PartialEq)]
pub struct RoundReport {
    /// The optimizer's reports, one per restart (a restart follows each remeshing that changed
    /// the mesh during the round).
    pub segments: Vec<OptimizeReport>,
    /// The segments merged: counters summed, `energy_initial` of the first and `energy_final` of
    /// the last segment, histories concatenated (a restart's first entry is the energy right
    /// after the remeshing), `mean_edge_length` of the first, `max_vertex_displacement` the
    /// largest over segments, and the status of the last segment (`MaxIterations` if the last
    /// segment was only stopped to remesh at the end of the budget).
    pub optimize: OptimizeReport,
    /// Every remeshing run in this round, in order: the one before the optimization (if enabled),
    /// then those triggered by the check. `energy_before/after` are filled in.
    pub remesh: Vec<RemeshReport>,
    /// Vertices at the end of the round.
    pub n_vertices: usize,
    /// Faces at the end of the round.
    pub n_faces: usize,
    /// Energy at the start of the round, after the subdivision (if any) and before any
    /// remeshing. Equal to the previous round's final energy when there is no subdivision (up
    /// to the energy of the new round's `EnergyOptions`).
    pub energy_after_subdivision: f64,
    /// Seconds spent subdividing.
    pub seconds_subdivide: f64,
    /// Seconds spent remeshing (including the energy evaluations for the report).
    pub seconds_remesh: f64,
    /// Seconds in the optimizer.
    pub seconds_optimize: f64,
}

/// Why a flow ended before a round (D32).
#[derive(Debug, Clone, PartialEq)]
pub struct FlowStop {
    /// The round that was not run.
    pub round: usize,
    /// Whether the undefined state is the result of that round's subdivision.
    pub after_subdivision: bool,
    /// The vertex where the energy is undefined (id in the state that would have started the round).
    pub vertex: u32,
    /// Why it is undefined (SPEC §3.6).
    pub reason: &'static str,
}

/// Report of [`flow`].
#[derive(Debug, Clone, PartialEq)]
pub struct FlowReport {
    /// One entry per round that ran.
    pub rounds: Vec<RoundReport>,
    /// Whether the progress callback ended the flow early (remaining rounds are not run).
    pub cancelled: bool,
    /// `Some` if the flow ended before a round because the energy of that round is undefined at
    /// its starting state (D32); the returned mesh is the state after the previous round.
    pub stopped: Option<FlowStop>,
    /// Smallest corner angle (radians) of the returned mesh. At least
    /// `optimize.min_face_angle` unless the input already was below it (then the flow failed).
    pub min_angle: f64,
    /// Faces of the returned mesh with a corner angle below the small angle of the remeshing
    /// (default 1 degree when remeshing is off): the ones the remeshing could not fix and the
    /// optimizer left (D31). The mesh is valid and admissible, but these faces are thin.
    pub bad_faces: usize,
    /// Wall time.
    pub seconds: f64,
}

/// Result of [`flow`].
#[derive(Debug, Clone, PartialEq)]
pub struct FlowResult {
    /// The final mesh.
    pub mesh: TriMesh,
    /// Its positions.
    pub positions: Vec<[f64; 3]>,
    /// For each final vertex, the id of the input vertex it descends from, or
    /// [`ddg_core::topology::NONE`] for a vertex created by subdivision. After a collapse
    /// the surviving vertex keeps its identity.
    pub vertex_source: Vec<u32>,
    /// For each final face, the input face it descends from. Subdivision: the four children of a
    /// face have its parent. Remeshing: a collapse removes faces and the others keep their
    /// parent; a flip rewrites two faces in place, the new face `(a, q, p)` keeps the parent of
    /// the flagged cap triangle and `(q, b, p)` that of its neighbour (SPEC §5.3).
    pub face_parent: Vec<u32>,
    /// The pinned vertices in the final numbering (ascending).
    pub pinned: Vec<u32>,
    /// What happened.
    pub report: FlowReport,
}

/// Number of faces with a corner angle below `small`.
fn count_bad(mesh: &TriMesh, pos: &[[f64; 3]], small: f64) -> usize {
    mesh.faces()
        .iter()
        .filter(|f| {
            let a = corner_angles(f.map(|v| pos[v as usize]));
            a[0].min(a[1]).min(a[2]) < small
        })
        .count()
}

fn merge_reports(segments: &[OptimizeReport], final_status: OptimizeStatus) -> OptimizeReport {
    let first = &segments[0];
    let last = &segments[segments.len() - 1];
    let mut m = OptimizeReport {
        status: final_status,
        iterations: 0,
        evaluations: 0,
        line_search_trials: 0,
        undefined_trials: 0,
        inadmissible_trials: 0,
        capped_steps: 0,
        weak_steps: 0,
        skipped_pairs: 0,
        memory_resets: 0,
        energy_initial: first.energy_initial,
        energy_final: last.energy_final,
        energy_history: Vec::new(),
        gradient_norm_history: Vec::new(),
        mean_edge_length: first.mean_edge_length,
        max_vertex_displacement: 0.0,
        energy_increases: 0,
        near_ties: last.near_ties,
        seconds: 0.0,
    };
    for s in segments {
        m.iterations += s.iterations;
        m.evaluations += s.evaluations;
        m.line_search_trials += s.line_search_trials;
        m.undefined_trials += s.undefined_trials;
        m.inadmissible_trials += s.inadmissible_trials;
        m.capped_steps += s.capped_steps;
        m.weak_steps += s.weak_steps;
        m.skipped_pairs += s.skipped_pairs;
        m.memory_resets += s.memory_resets;
        m.energy_history.extend_from_slice(&s.energy_history);
        m.gradient_norm_history.extend_from_slice(&s.gradient_norm_history);
        m.max_vertex_displacement = m.max_vertex_displacement.max(s.max_vertex_displacement);
        m.energy_increases += s.energy_increases;
        m.seconds += s.seconds;
    }
    m
}

struct State {
    mesh: TriMesh,
    pos: Vec<[f64; 3]>,
    source: Vec<u32>,
    face_parent: Vec<u32>,
    pinned: Vec<u32>,
}

/// Runs the developability flow (SPEC §5.1).
///
/// The `progress` callback is called after every accepted iteration with the round index and the
/// iteration count of the round (over restarts); its `positions` belong to the mesh of the
/// current moment (the topology can change between calls, when `vertex_source` would be needed to
/// follow vertices). Returning [`Control::Break`] ends the flow with the last accepted state; the
/// report says `cancelled`.
///
/// # Errors
/// Invalid options, wrong lengths, non-finite positions, pinned ids out of range, a degenerate
/// face or an undefined energy at a starting state, and anything [`optimize`], [`subdivide`] and
/// [`remesh`] report.
pub fn flow(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &FlowOptions,
    progress: Option<&mut dyn FnMut(&Progress) -> Control>,
) -> Result<FlowResult, DevError> {
    opts.validate()?;
    let nv = mesh.n_vertices();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch { what: "positions", expected: nv, got: positions.len() });
    }
    if let Some(v) = positions.iter().position(|p| !v3::is_finite(*p)) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    if let Some(&p) = opts.pinned.iter().find(|&&p| p as usize >= nv) {
        return Err(DevError::VertexOutOfRange { what: "pinned", vertex: p, n_vertices: nv });
    }
    match opts.threads {
        None => flow_impl(mesh, positions, opts, progress, false),
        Some(n) => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| DevError::Internal(format!("cannot build a thread pool of {n}: {e}")))?;
            let Some(user) = progress else {
                return pool.install(|| flow_impl(mesh, positions, opts, None, true));
            };
            // The flow runs in the pool; the (not necessarily `Send`) callback stays on this
            // thread and is served over channels: every progress report is copied and answered.
            let (req_tx, req_rx) = mpsc::channel::<Msg>();
            let (resp_tx, resp_rx) = mpsc::channel::<Control>();
            let mut result = None;
            std::thread::scope(|s| {
                let done_tx = req_tx.clone();
                let pool = &pool;
                s.spawn(move || {
                    let r = pool.install(move || {
                        let mut forward = |p: &Progress| -> Control {
                            let m = Msg::Progress {
                                round: p.round,
                                iteration: p.iteration,
                                energy: p.energy,
                                gradient_norm: p.gradient_norm,
                                positions: p.positions.to_vec(),
                            };
                            if req_tx.send(m).is_err() {
                                return Control::Break;
                            }
                            resp_rx.recv().unwrap_or(Control::Break)
                        };
                        flow_impl(mesh, positions, opts, Some(&mut forward), true)
                    });
                    let _ = done_tx.send(Msg::Done(Box::new(r)));
                });
                while let Ok(m) = req_rx.recv() {
                    match m {
                        Msg::Progress { round, iteration, energy, gradient_norm, positions } => {
                            let c = user(&Progress { round, iteration, energy, gradient_norm, positions: &positions });
                            let _ = resp_tx.send(c);
                        }
                        Msg::Done(r) => {
                            result = Some(*r);
                            break;
                        }
                    }
                }
            });
            result.unwrap_or_else(|| Err(DevError::Internal("the flow thread ended without a result".into())))
        }
    }
}

enum Msg {
    Progress { round: usize, iteration: usize, energy: f64, gradient_norm: f64, positions: Vec<[f64; 3]> },
    Done(Box<Result<FlowResult, DevError>>),
}

fn flow_impl(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &FlowOptions,
    mut progress: Option<&mut dyn FnMut(&Progress) -> Control>,
    in_pool: bool,
) -> Result<FlowResult, DevError> {
    let start = Instant::now();
    let mut pinned: Vec<u32> = opts.pinned.clone();
    pinned.sort_unstable();
    pinned.dedup();
    let mut st = State {
        mesh: mesh.clone(),
        pos: positions.to_vec(),
        source: (0..mesh.n_vertices() as u32).collect(),
        face_parent: (0..mesh.n_faces() as u32).collect(),
        pinned,
    };
    let mut rounds: Vec<RoundReport> = Vec::new();
    let mut cancelled = false;
    let mut stopped: Option<FlowStop> = None;

    for (ri, round) in opts.rounds.iter().enumerate() {
        let mut energy = round.energy.clone();
        if in_pool {
            energy.threads = None; // use the pool we are in
        }
        let mut seconds_remesh = 0.0;

        // 1. Subdivision. The energy of the round must be defined at its starting state (the
        // intrinsic energies are not on the midpoint of a flat fold, D32): if it is not, the flow
        // stops before the round, with the state of the previous round, and says so (`stopped`).
        let t = Instant::now();
        let candidate = if round.subdivide_before {
            let sd = subdivide(&st.mesh, &st.pos)?;
            let pinned = sd.pinned_after(&st.mesh, &st.pinned)?;
            Some((sd, pinned))
        } else {
            None
        };
        let start_state = candidate.as_ref().map_or((&st.mesh, &st.pos[..]), |(sd, _)| (&sd.mesh, &sd.positions[..]));
        let energy_start = match evaluate_energy(start_state.0, start_state.1, &energy) {
            Ok((e, _)) => e,
            Err(DevError::UndefinedEnergy { vertex, reason }) if ri > 0 => {
                stopped = Some(FlowStop { round: ri, after_subdivision: round.subdivide_before, vertex, reason });
                break;
            }
            Err(e) => return Err(e),
        };
        if let Some((sd, pinned)) = candidate {
            st.pinned = pinned;
            st.source.resize(sd.mesh.n_vertices(), NONE);
            st.face_parent = sd.face_parent.iter().map(|&p| st.face_parent[p as usize]).collect();
            st.mesh = sd.mesh;
            st.pos = sd.positions;
        }
        let seconds_subdivide = if round.subdivide_before { t.elapsed().as_secs_f64() } else { 0.0 };

        // 2. Remesh before optimizing.
        let mut remesh_reports: Vec<RemeshReport> = Vec::new();
        let mut known_unfixable = 0usize;
        if let Some(ro) = &opts.remesh {
            let t = Instant::now();
            let rr = run_remesh(&mut st, ro, &energy, opts.optimize.min_face_angle)?;
            remesh_reports.push(rr);
            known_unfixable = count_bad(&st.mesh, &st.pos, ro.small_angle);
            seconds_remesh += t.elapsed().as_secs_f64();
        }

        // 3. Optimize with the remeshing check.
        let mut oo = opts.optimize.clone();
        let mut remaining = round.max_iterations;
        let mut done = 0usize;
        let mut segments: Vec<OptimizeReport> = Vec::new();
        let mut seconds_optimize = 0.0;
        let final_status: OptimizeStatus;
        let mut zero_runs = 0usize;
        loop {
            oo.max_iterations = remaining;
            oo.pinned = st.pinned.clone();
            let mut need_remesh = false;
            let mut user_break = false;
            let small = opts.remesh.as_ref().map(|r| r.small_angle);
            let interval = opts.remesh.as_ref().map_or(1, |r| r.check_interval);
            let base = done;
            let mesh_ref = &st.mesh;
            let res = {
                let mut cb = |p: &Progress| -> Control {
                    if let Some(user) = progress.as_deref_mut() {
                        let q = Progress { round: ri, iteration: base + p.iteration, ..*p };
                        if user(&q) == Control::Break {
                            user_break = true;
                            return Control::Break;
                        }
                    }
                    if let Some(small) = small
                        && p.iteration.is_multiple_of(interval)
                        && count_bad(mesh_ref, p.positions, small) > known_unfixable
                    {
                        need_remesh = true;
                        return Control::Break;
                    }
                    Control::Continue
                };
                let t = Instant::now();
                let r = optimize(&st.mesh, &st.pos, &energy, &oo, Some(&mut cb));
                seconds_optimize += t.elapsed().as_secs_f64();
                r?
            };
            st.pos = res.positions;
            let rep = res.report;
            let iters = rep.iterations;
            let status = rep.status;
            segments.push(rep);
            done += iters;
            remaining = remaining.saturating_sub(iters);
            if user_break {
                cancelled = true;
                final_status = OptimizeStatus::Cancelled;
                break;
            }
            let failed = matches!(status, OptimizeStatus::LineSearchFailed | OptimizeStatus::Stationary | OptimizeStatus::NonDescent);
            if !(need_remesh || (failed && opts.remesh.is_some())) {
                final_status = status;
                break;
            }
            // remesh (a failure may be caused by a bad triangle the check has not flagged)
            let ro = opts.remesh.as_ref().expect("remesh requested only when enabled");
            let t = Instant::now();
            let rr = run_remesh(&mut st, ro, &energy, opts.optimize.min_face_angle)?;
            seconds_remesh += t.elapsed().as_secs_f64();
            let changed = rr.flips + rr.collapses > 0 && rr.rejected.is_none();
            remesh_reports.push(rr);
            known_unfixable = count_bad(&st.mesh, &st.pos, ro.small_angle);
            zero_runs = if iters == 0 { zero_runs + 1 } else { 0 };
            if remaining == 0 {
                final_status = if need_remesh { OptimizeStatus::MaxIterations } else { status };
                break;
            }
            if zero_runs > 100 || (!changed && !need_remesh) {
                final_status = status;
                break;
            }
            // changed: restart on the new mesh; flagged but nothing fixable: resume with the
            // current bad triangles as the baseline (each restart makes progress, so this ends)
        }
        let optimize_report = merge_reports(&segments, final_status);
        rounds.push(RoundReport {
            segments,
            optimize: optimize_report,
            remesh: remesh_reports,
            n_vertices: st.mesh.n_vertices(),
            n_faces: st.mesh.n_faces(),
            energy_after_subdivision: energy_start,
            seconds_subdivide,
            seconds_remesh,
            seconds_optimize,
        });
        if cancelled {
            break;
        }
    }
    let small = opts.remesh.as_ref().map_or(RemeshOptions::default().small_angle, |r| r.small_angle);
    let min_angle = st
        .mesh
        .faces()
        .iter()
        .flat_map(|f| corner_angles(f.map(|v| st.pos[v as usize])))
        .fold(f64::INFINITY, f64::min);
    let bad_faces = count_bad(&st.mesh, &st.pos, small);
    Ok(FlowResult {
        mesh: st.mesh,
        positions: st.pos,
        vertex_source: st.source,
        face_parent: st.face_parent,
        pinned: st.pinned,
        report: FlowReport { rounds, cancelled, stopped, min_angle, bad_faces, seconds: start.elapsed().as_secs_f64() },
    })
}

/// Remeshes the state in place, keeping the id bookkeeping, and fills in the report's energies.
/// A result that the optimizer could not start from (undefined energy, a face below `floor`) is
/// not applied (`rejected`, D32); everything else is unchanged by that.
fn run_remesh(
    st: &mut State,
    ro: &RemeshOptions,
    energy: &EnergyOptions,
    floor: f64,
) -> Result<RemeshReport, DevError> {
    let before = evaluate_energy(&st.mesh, &st.pos, energy).ok().map(|e| e.0);
    let rm = remesh(&st.mesh, &st.pos, ro, &st.pinned)?;
    let mut report = rm.report.clone();
    report.energy_before = before;
    if rm.changed() {
        let after = evaluate_energy(&rm.mesh, &rm.positions, energy);
        let reason = match &after {
            Err(DevError::UndefinedEnergy { .. }) => Some("the energy is undefined on the result"),
            Err(_) => Some("the energy cannot be evaluated on the result"),
            Ok(_) => check_admissible(&rm.mesh, &rm.positions, floor, START_SLACK)
                .err()
                .map(|_| "a face of the result is below the angle floor"),
        };
        if let Some(r) = reason {
            report.rejected = Some(r);
            report.energy_after = before;
            return Ok(report);
        }
        st.source = rm.vertex_to_old.iter().map(|&o| st.source[o as usize]).collect();
        st.face_parent = rm.face_to_old.iter().map(|&o| st.face_parent[o as usize]).collect();
        st.pinned = rm.pinned;
        st.mesh = rm.mesh;
        st.pos = rm.positions;
        report.energy_after = after.ok().map(|e| e.0);
    } else {
        report.energy_after = before;
    }
    Ok(report)
}
