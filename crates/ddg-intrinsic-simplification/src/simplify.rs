//! Intrinsic coarsening (paper Sec. 6, Alg. 1; SPEC §5): the priority queue of ICE costs,
//! targets, fixed vertices, snapshots and the report; and the option set that reproduces the
//! reference's choices ([`ReferenceCompat`], SPEC §12.2).

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BinaryHeap};
use std::time::Instant;

use ddg_intrinsic::{DelaunayCriterion, DelaunayStatus};

use crate::error::{IemError, RemovalRefusal};
use crate::flip::FlippedEdge;
use crate::mesh::{CoarseMesh, PointLocation};
use crate::metric::{ChannelState, ChannelWeights, Frames, MetricOptions};
use crate::removal::{PreparedRemoval, RemovalKind};

/// When [`simplify`] stops (SPEC §5 item 4). It also stops when nothing is removable (the
/// minimum of the queue is `∞`, S8) or the queue is empty.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Target {
    /// At this many vertices in total (the paper's `n`).
    Vertices(usize),
    /// At `round(r · n₀)` vertices, `n₀` the vertex count when `simplify` starts, `r ∈ [0, 1]`.
    Ratio(f64),
    /// Each connected component at `round(r · n_c)` vertices, `n_c` its vertex count when
    /// `simplify` starts (paper Sec. 8.2 normalizes by components), `r ∈ [0, 1]`.
    PerComponent(f64),
    /// When the cheapest removal costs more than this (`≥ 0`).
    MaxCost(f64),
    /// Only when nothing is removable.
    AllRemovable,
}

/// Which vertices get a new cost after a removal (S6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Requeue {
    /// `N_i` plus the four vertices of every flip of the removal (degree, fallback, ear and
    /// Delaunay flips): every vertex whose star or masses changed (S6, default).
    StarsChanged,
    /// `N_i` only (paper Alg. 1 lines 13–15; the reference).
    OneRing,
}

/// The reference's choices that are not bugs (SPEC §12.2), for differential tests against
/// the reference driver; see [`SimplifyOptions::reference_compat`]. Its fields cover what
/// [`SimplifyOptions`] and [`MetricOptions`] do not.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceCompat {
    /// Multiply the channel weights into the masses (R2).
    pub fold_weights: bool,
    /// The frames of the queue key (R16: before flattening).
    pub key_frames: Frames,
    /// The frames of the commit (R16, R19: before flattening at `i`, after at `j`).
    pub commit_frames: Frames,
    /// Seed the Delaunay flips after a removal with every edge of the faces within two edge
    /// steps of the merged face (R21) instead of the edges among `N_i` (S15).
    pub delaunay_two_ring: bool,
    /// The Delaunay test's tolerance on the opposite-angle sum (R22: `1e-5`); `None` keeps the
    /// mesh's criterion.
    pub delaunay_tolerance: Option<f64>,
    /// The flattening's relative tolerance (the reference stops at `|K| < 1e-7`, i.e.
    /// `1e-7/2π` relative to `2π`); `None` keeps the mesh's.
    pub flatten_tolerance: Option<f64>,
    /// Break ties of equal costs by the removal count at which the entry was pushed, then by
    /// vertex (the reference's `(cost, time stamp, v)`), instead of by vertex.
    pub stamp_ties: bool,
}

impl Default for ReferenceCompat {
    fn default() -> Self {
        ReferenceCompat {
            fold_weights: true,
            key_frames: Frames::BeforeFlattening,
            commit_frames: Frames::ReferenceCommit,
            delaunay_two_ring: true,
            delaunay_tolerance: Some(1e-5),
            flatten_tolerance: Some(1e-7 / std::f64::consts::TAU),
            stamp_ties: true,
        }
    }
}

/// The options of [`simplify`] (SPEC §5).
#[derive(Debug, Clone, PartialEq)]
pub struct SimplifyOptions {
    /// When to stop.
    pub target: Target,
    /// The channels and their weights (S10: the curvature channels, weight 1).
    pub channels: ChannelWeights,
    /// The cost (memory-based by default).
    pub metric: MetricOptions,
    /// Flip to intrinsic Delaunay before computing masses and costs (S3, Alg. 1 line 3;
    /// default on; the reference skips it, R1).
    pub initial_delaunay: bool,
    /// Which vertices get a new cost after a removal (S6).
    pub requeue: Requeue,
    /// Fine vertices to fix (paper Sec. 8.6.3), in addition to those the mesh already has
    /// fixed. They stay fixed on the mesh afterwards.
    pub fixed: Vec<u32>,
    /// Vertex counts (strictly descending) at which to record a [`Snapshot`] (S24).
    pub snapshots: Vec<usize>,
    /// Record the `(vertex, key, cost)` of every removal in the report.
    pub record_sequence: bool,
    /// Check that a queued key equals the commit cost whenever nothing it depends on changed
    /// since it was computed (S5), counted in [`KeyCheck`]. Cheap; on by default.
    pub check_keys: bool,
    /// The reference's choices (SPEC §12.2); `None` (default) for the crate's.
    pub reference: Option<ReferenceCompat>,
}

impl Default for SimplifyOptions {
    fn default() -> Self {
        SimplifyOptions {
            target: Target::AllRemovable,
            channels: ChannelWeights::curvature(),
            metric: MetricOptions::default(),
            initial_delaunay: true,
            requeue: Requeue::StarsChanged,
            fixed: Vec::new(),
            snapshots: Vec::new(),
            record_sequence: false,
            check_keys: true,
            reference: None,
        }
    }
}

impl SimplifyOptions {
    /// The options with which this crate reproduces the reference's non-bug choices (SPEC
    /// §12.2): no initial Delaunay (R1), `mixture(w)` with the weights folded into the masses
    /// (R2), the reference's boundary normalization `2π/(Θ + π)` in the transport (R3), the cost
    /// over receiving neighbours only (R11), equal shares only for an exactly zero change with
    /// `+1e-6` (R12), keys in the frames before flattening and the commit in the frames before
    /// flattening at `i` (R16, R19), `N_i` re-queued (Alg. 1), 2-ring Delaunay seeding (R21),
    /// the Delaunay tolerance `1e-5` (R22), the flattening tolerance `|K| < 1e-7`, and ties by
    /// time stamp. The boundary removal by degree 2 and a merge gives the same mesh as the
    /// reference's flip to degree 1 and ear deletion (R9; see SPEC §12.2), so there is no
    /// option for it.
    pub fn reference_compat(w: f64, target: Target) -> SimplifyOptions {
        SimplifyOptions {
            target,
            channels: ChannelWeights::mixture(w),
            metric: MetricOptions {
                memory: true,
                flat: crate::metric::FlatRule::ReferenceEpsilon,
                receiving_only: true,
                boundary_frames: crate::metric::BoundaryFrames::Reference,
            },
            initial_delaunay: false,
            requeue: Requeue::OneRing,
            reference: Some(ReferenceCompat::default()),
            ..SimplifyOptions::default()
        }
    }
}

/// Why [`simplify`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    /// The target vertex count (total or every component) was reached.
    TargetReached,
    /// The cheapest removal costs more than [`Target::MaxCost`].
    MaxCost,
    /// The cheapest removal costs `∞`: nothing is removable (S8).
    NothingRemovable,
    /// The queue ran empty.
    QueueEmpty,
}

/// One removal of a [`simplify`] run (for differential tests).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemovalRecord {
    /// The removed vertex (fine id).
    pub vertex: u32,
    /// Its kind.
    pub kind: RemovalKind,
    /// The key it was popped with.
    pub key: f64,
    /// The cost of the committed removal.
    pub cost: f64,
}

/// The check of S5: a key equals the commit cost when nothing it depends on changed.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct KeyCheck {
    /// Committed removals whose key was computed after the last change of the vertex and of
    /// all its neighbours (star or masses).
    pub n_clean: usize,
    /// Among those, removals whose key differs from the commit cost (bitwise). Always 0 unless
    /// the key and commit frames differ ([`ReferenceCompat`], R16).
    pub n_clean_mismatch: usize,
    /// Removals whose key differs from the commit cost (a stale key: the masses of a 2-ring
    /// vertex changed, or the reference's frames, R16).
    pub n_differ: usize,
    /// The largest `|key − cost| / max(|cost|, tiny)` over the removals.
    pub max_relative_difference: f64,
}

/// Seconds per phase of [`simplify`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SimplifyTimings {
    /// The initial Delaunay flips.
    pub initial_delaunay: f64,
    /// The initial masses.
    pub channels: f64,
    /// The initial cost of every vertex.
    pub initial_costs: f64,
    /// The committed removals (flattening, flips, removal, points, transport, local Delaunay).
    pub removals: f64,
    /// The new costs after removals.
    pub requeue: f64,
    /// The snapshots.
    pub snapshots: f64,
    /// Everything.
    pub total: f64,
}

/// A connected component's vertex counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComponentCount {
    /// Vertices when `simplify` started.
    pub before: usize,
    /// Vertices at the end.
    pub after: usize,
}

/// The coarse mesh at a given vertex count during [`simplify`] (S24), compacted: vertex `k` is
/// fine vertex `vertices[k]`; faces and lengths in the coarse mesh's face order. With the
/// locations of the previous level's vertices it gives the prolongation between levels.
#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    /// The removals done when it was taken.
    pub n_removed: usize,
    /// Its vertices as fine ids.
    pub vertices: Vec<u32>,
    /// Its faces (indices into `vertices`; a Δ-face may repeat one).
    pub faces: Vec<[u32; 3]>,
    /// Per face, `l[k]` the length of the side from slot `k` to slot `k + 1`.
    pub lengths: Vec<[f64; 3]>,
    /// The vertices of the previous level (the previous snapshot; the vertices when `simplify`
    /// started for the first), as fine ids.
    pub previous_vertices: Vec<u32>,
    /// Where each of them is in this snapshot (faces are indices into `faces`, corners fine
    /// ids).
    pub previous: Vec<PointLocation>,
}

/// What [`simplify`] did (SPEC §5).
#[derive(Debug, Clone, PartialEq)]
pub struct SimplifyReport {
    /// Why it stopped.
    pub stop: StopReason,
    /// Vertices before and after.
    pub n_vertices_before: usize,
    /// Vertices after.
    pub n_vertices_after: usize,
    /// The initial Delaunay flips (S3), if any.
    pub n_initial_delaunay_flips: usize,
    /// How the initial Delaunay pass ended.
    pub initial_delaunay_status: Option<DelaunayStatus>,
    /// Removals of interior vertices.
    pub n_interior: usize,
    /// Removals of regular boundary vertices.
    pub n_boundary: usize,
    /// Removals of ears.
    pub n_ear: usize,
    /// The vertices left at the end whose last evaluation was `∞`, by reason (the skips).
    pub skips: BTreeMap<&'static str, usize>,
    /// Every evaluation (initial and re-queue) that gave `∞`, by reason.
    pub infinite_evaluations: BTreeMap<&'static str, usize>,
    /// Popped vertices with a finite key whose removal failed at commit (reverted, not
    /// re-queued until a neighbour changes; S7), by reason.
    pub commit_failures: BTreeMap<&'static str, usize>,
    /// Entries popped from the queue.
    pub n_pops: usize,
    /// Popped entries that were out of date (lazy deletion).
    pub n_stale_pops: usize,
    /// Popped entries of components that had reached their target.
    pub n_component_done_pops: usize,
    /// Tentative removals evaluated.
    pub n_evaluations: usize,
    /// The largest number of flattening steps of a committed removal.
    pub flatten_iterations_max: usize,
    /// The mean number of flattening steps over committed removals.
    pub flatten_iterations_mean: f64,
    /// Springborn's fallback flips in committed removals.
    pub n_fallback_flips: usize,
    /// Flips to degree three or two in committed removals (with the ear flips).
    pub n_degree_flips: usize,
    /// Delaunay flips after removals.
    pub n_delaunay_flips: usize,
    /// Delaunay flips refused (non-Delaunay edges that could not be flipped).
    pub n_delaunay_refused: usize,
    /// Points clamped onto a side during the run (S16).
    pub n_clamped: u64,
    /// Per connected component (in order of the smallest fine id), its vertex counts.
    pub components: Vec<ComponentCount>,
    /// The check of S5.
    pub keys: KeyCheck,
    /// Seconds per phase.
    pub timings: SimplifyTimings,
    /// Every removal, if asked ([`SimplifyOptions::record_sequence`]).
    pub sequence: Vec<RemovalRecord>,
    /// The snapshots reached, in order.
    pub snapshots: Vec<Snapshot>,
    /// The masses and vectors at the end.
    pub channels: ChannelState,
}

/// A queue entry: ordered by cost, then (with stamps) push time, then vertex.
#[derive(Debug, Clone, Copy)]
struct Entry {
    cost: f64,
    stamp: u32,
    fine: u32,
    version: u32,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Entry {}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cost
            .total_cmp(&other.cost)
            .then(self.stamp.cmp(&other.stamp))
            .then(self.fine.cmp(&other.fine))
            .then(self.version.cmp(&other.version))
    }
}

/// The name of a refusal reason in the report's counts.
pub fn refusal_name(r: &RemovalRefusal) -> &'static str {
    match r {
        RemovalRefusal::NotAVertex => "not_a_vertex",
        RemovalRefusal::Fixed => "fixed",
        RemovalRefusal::SelfFace { .. } => "self_face",
        RemovalRefusal::BoundarySelfEdge => "boundary_self_edge",
        RemovalRefusal::InteriorDegree { .. } => "interior_degree",
        RemovalRefusal::EarNotFlippable { .. } => "ear_not_flippable",
        RemovalRefusal::FlattenFailed { .. } => "flatten_failed",
        RemovalRefusal::NotFlattenable { .. } => "not_flattenable",
        RemovalRefusal::Stuck { .. } => "stuck",
        RemovalRefusal::SelfEdge => "self_edge",
        RemovalRefusal::Degenerate { .. } => "degenerate",
        RemovalRefusal::Internal { .. } => "internal",
    }
}

/// The state of one run.
struct Run<'a> {
    mesh: &'a mut CoarseMesh,
    options: &'a SimplifyOptions,
    state: ChannelState,
    key_frames: Frames,
    commit_frames: Frames,
    stamp_ties: bool,
    heap: BinaryHeap<Reverse<Entry>>,
    version: Vec<u32>,
    /// Per fine id: the last refusal of its evaluation (`None` if finite).
    last_refusal: Vec<Option<&'static str>>,
    /// Per fine id: the removal count at its last change (star or masses), and when its key
    /// was computed (S5 check).
    touched: Vec<u32>,
    key_step: Vec<u32>,
    n_removed: u32,
    alphas: Vec<f64>,
    report: SimplifyReport,
}

impl Run<'_> {
    /// Evaluates the cost of removing `v` now and queues it.
    fn push(&mut self, v: u32) {
        self.report.n_evaluations += 1;
        let cost = match self.mesh.tentative_removal(v) {
            Ok(prepared) => {
                self.last_refusal[v as usize] = None;
                let c = self.state.evaluate(
                    &prepared,
                    &self.options.metric,
                    self.key_frames,
                    &mut self.alphas,
                    None,
                );
                if c.is_nan() { f64::INFINITY } else { c }
            }
            Err(f) => {
                let name = refusal_name(&f.reason);
                self.last_refusal[v as usize] = Some(name);
                *self.report.infinite_evaluations.entry(name).or_default() += 1;
                f64::INFINITY
            }
        };
        let version = self.version[v as usize].wrapping_add(1);
        self.version[v as usize] = version;
        self.key_step[v as usize] = self.n_removed;
        self.heap.push(Reverse(Entry {
            cost,
            stamp: if self.stamp_ties { self.n_removed } else { 0 },
            fine: v,
            version,
        }));
    }

    fn touch(&mut self, v: u32) {
        if let Some(t) = self.touched.get_mut(v as usize) {
            *t = self.n_removed;
        }
    }
}

/// The faces within two edge steps of face `f` (R21).
fn two_ring_edges(mesh: &CoarseMesh, f: usize) -> Vec<u32> {
    let c = mesh.complex();
    let mut faces = vec![f];
    let mut start = 0;
    for _ in 0..2 {
        let end = faces.len();
        for k in start..end {
            let g = faces[k];
            for s in 0..3 {
                if let Some(t) = c.twin(3 * g + s) {
                    let h = t / 3;
                    if !faces.contains(&h) {
                        faces.push(h);
                    }
                }
            }
        }
        start = end;
    }
    let mut edges: Vec<u32> = faces
        .iter()
        .flat_map(|&g| (0..3).map(move |s| c.edge(3 * g + s) as u32))
        .collect();
    edges.sort_unstable();
    edges.dedup();
    edges
}

/// The compacted coarse mesh with the locations of the previous level's vertices.
fn take_snapshot(mesh: &CoarseMesh, n_removed: usize, previous_vertices: Vec<u32>) -> Snapshot {
    let c = mesh.complex();
    let lengths = mesh.edge_lengths();
    let previous = previous_vertices
        .iter()
        .map(|&v| {
            mesh.point_location(v as usize)
                .unwrap_or(PointLocation::Vertex(v))
        })
        .collect();
    Snapshot {
        n_removed,
        vertices: mesh.vertices().to_vec(),
        faces: c.faces().to_vec(),
        lengths: (0..c.n_faces())
            .map(|f| [0, 1, 2].map(|k| lengths[c.edge(3 * f + k)]))
            .collect(),
        previous_vertices,
        previous,
    }
}

fn check_options(mesh: &CoarseMesh, o: &SimplifyOptions) -> Result<(), IemError> {
    let bad = |what: String| Err(IemError::InvalidOptions { what });
    match o.target {
        Target::Ratio(r) | Target::PerComponent(r) if !(0.0..=1.0).contains(&r) => {
            return bad(format!("the target ratio {r} is not in [0, 1]"));
        }
        Target::MaxCost(c) if c.is_nan() || c < 0.0 => {
            return bad(format!("the maximum cost {c} is not a non-negative number"));
        }
        _ => {}
    }
    if let Some(w) = o.snapshots.windows(2).find(|w| w[0] <= w[1]) {
        return bad(format!(
            "snapshot counts must be strictly descending, but {} is followed by {}",
            w[0], w[1]
        ));
    }
    if let Some(&v) = o
        .fixed
        .iter()
        .find(|&&v| v as usize >= mesh.n_fine_vertices())
    {
        return bad(format!(
            "fixed vertex {v} is out of range (there are {} fine vertices)",
            mesh.n_fine_vertices()
        ));
    }
    if let Some(r) = &o.reference {
        if let Some(t) = r.delaunay_tolerance
            && !(t.is_finite() && t >= 0.0)
        {
            return bad(format!("the Delaunay tolerance {t} is not finite and >= 0"));
        }
        if let Some(t) = r.flatten_tolerance
            && !(t.is_finite() && t > 0.0)
        {
            return bad(format!(
                "the flattening tolerance {t} is not finite and > 0"
            ));
        }
    }
    Ok(())
}

/// The connected components of the coarse mesh's vertices: per fine id its component (`u32::MAX`
/// for removed vertices), numbered in order of the smallest fine id.
fn components(mesh: &CoarseMesh) -> (Vec<u32>, usize) {
    let c = mesh.complex();
    let n = c.n_vertices();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for f in c.faces() {
        for k in 0..3 {
            let (a, b) = (
                find(&mut parent, f[k] as usize),
                find(&mut parent, f[(k + 1) % 3] as usize),
            );
            if a != b {
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    // Order components by the smallest fine id in them.
    let mut order: Vec<(u32, usize)> = (0..n).map(|v| (mesh.fine_vertex(v), v)).collect();
    order.sort_unstable();
    let mut label = vec![u32::MAX; n];
    let mut n_comp = 0;
    let mut out = vec![u32::MAX; mesh.fine_id_bound()];
    for (fine, v) in order {
        let r = find(&mut parent, v);
        if label[r] == u32::MAX {
            label[r] = n_comp as u32;
            n_comp += 1;
        }
        out[fine as usize] = label[r];
    }
    (out, n_comp)
}

/// Coarsens `mesh` by the ICE metric (paper Alg. 1; SPEC §5): initial Delaunay flips (S3),
/// masses from the channels and a tentative cost per vertex in a queue keyed `(cost, fine id)`
/// with lazy deletion; then the cheapest vertex is removed, its mass transported to its
/// neighbours, the mesh flipped back to Delaunay locally, and the vertices whose star or masses
/// changed get new costs (S6), until the target or until the cheapest cost is `∞` (S8). A
/// removal that fails at commit is reverted and the vertex is not re-queued until a neighbour
/// changes (S7).
///
/// # Errors
/// [`IemError::RemovalPending`]; [`IemError::InvalidOptions`] (see the options); and
/// [`IemError::Internal`] from a flip or removal that meets a state its invariants exclude
/// (never expected; the mesh may then be partly simplified).
pub fn simplify(
    mesh: &mut CoarseMesh,
    options: &SimplifyOptions,
) -> Result<SimplifyReport, IemError> {
    if mesh.removal_pending() {
        return Err(IemError::RemovalPending);
    }
    check_options(mesh, options)?;
    let saved_options = mesh.options;
    let compat = options.reference;
    if let Some(r) = compat {
        if let Some(t) = r.delaunay_tolerance {
            mesh.options.delaunay = DelaunayCriterion::AngleSum { tolerance: t };
        }
        if let Some(t) = r.flatten_tolerance {
            mesh.options.flatten_tolerance = t;
        }
        if r.delaunay_two_ring {
            mesh.options.delaunay_after_removal = false;
        }
    }
    let result = run(mesh, options);
    mesh.options = saved_options;
    result
}

fn run(mesh: &mut CoarseMesh, options: &SimplifyOptions) -> Result<SimplifyReport, IemError> {
    let t_start = Instant::now();
    let compat = options.reference;
    for &v in &options.fixed {
        mesh.set_fixed(v, true)?;
    }
    let n_before = mesh.n_vertices();
    // Per-vertex state is indexed by fine id, pre-refinement's inserted vertices included.
    let n_fine = mesh.fine_id_bound();
    let (component, n_components) = components(mesh);
    let mut comp_count = vec![0usize; n_components];
    for &v in mesh.vertices() {
        comp_count[component[v as usize] as usize] += 1;
    }
    let comp_before = comp_count.clone();

    // 1. Initial Delaunay (S3).
    let t0 = Instant::now();
    let initial = if options.initial_delaunay {
        Some(mesh.flip_to_delaunay()?)
    } else {
        None
    };
    let t_delaunay = t0.elapsed().as_secs_f64();

    // 2. Masses.
    let t0 = Instant::now();
    let fold = compat.is_some_and(|r| r.fold_weights);
    let state = ChannelState::new(mesh, &options.channels, fold)?;
    let t_channels = t0.elapsed().as_secs_f64();

    let report = SimplifyReport {
        stop: StopReason::QueueEmpty,
        n_vertices_before: n_before,
        n_vertices_after: n_before,
        n_initial_delaunay_flips: initial.as_ref().map_or(0, |r| r.flips.len()),
        initial_delaunay_status: initial.as_ref().map(|r| r.status),
        n_interior: 0,
        n_boundary: 0,
        n_ear: 0,
        skips: BTreeMap::new(),
        infinite_evaluations: BTreeMap::new(),
        commit_failures: BTreeMap::new(),
        n_pops: 0,
        n_stale_pops: 0,
        n_component_done_pops: 0,
        n_evaluations: 0,
        flatten_iterations_max: 0,
        flatten_iterations_mean: 0.0,
        n_fallback_flips: 0,
        n_degree_flips: 0,
        n_delaunay_flips: 0,
        n_delaunay_refused: 0,
        n_clamped: 0,
        components: Vec::new(),
        keys: KeyCheck::default(),
        timings: SimplifyTimings {
            initial_delaunay: t_delaunay,
            channels: t_channels,
            ..SimplifyTimings::default()
        },
        sequence: Vec::new(),
        snapshots: Vec::new(),
        channels: ChannelState::empty(),
    };
    let clamped_before = mesh.n_clamped_points();
    let mut run = Run {
        options,
        state,
        key_frames: compat.map_or(Frames::AfterFlattening, |r| r.key_frames),
        commit_frames: compat.map_or(Frames::AfterFlattening, |r| r.commit_frames),
        stamp_ties: compat.is_some_and(|r| r.stamp_ties),
        heap: BinaryHeap::with_capacity(4 * n_before),
        version: vec![0; n_fine],
        last_refusal: vec![None; n_fine],
        touched: vec![0; n_fine],
        key_step: vec![0; n_fine],
        n_removed: 0,
        alphas: Vec::with_capacity(16),
        report,
        mesh,
    };

    // Snapshots (S24).
    let mut next_snapshot = 0;
    let mut previous_vertices: Vec<u32> = run.mesh.vertices().to_vec();
    previous_vertices.sort_unstable();
    let mut take_snapshots = |run: &mut Run, previous_vertices: &mut Vec<u32>| {
        let t = Instant::now();
        while next_snapshot < options.snapshots.len()
            && run.mesh.n_vertices() <= options.snapshots[next_snapshot]
        {
            let snap = take_snapshot(
                run.mesh,
                run.n_removed as usize,
                std::mem::take(previous_vertices),
            );
            *previous_vertices = snap.vertices.clone();
            previous_vertices.sort_unstable();
            run.report.snapshots.push(snap);
            next_snapshot += 1;
        }
        run.report.timings.snapshots += t.elapsed().as_secs_f64();
    };
    take_snapshots(&mut run, &mut previous_vertices);

    // 3. Initial costs, in order of fine id.
    let t0 = Instant::now();
    let mut live: Vec<u32> = run.mesh.vertices().to_vec();
    live.sort_unstable();
    for &v in &live {
        run.push(v);
    }
    run.report.timings.initial_costs = t0.elapsed().as_secs_f64();

    // Targets.
    let total_target = match options.target {
        Target::Vertices(n) => n,
        Target::Ratio(r) => (r * n_before as f64).round() as usize,
        _ => 0,
    };
    let comp_target: Vec<usize> = match options.target {
        Target::PerComponent(r) => comp_before
            .iter()
            .map(|&n| (r * n as f64).round() as usize)
            .collect(),
        _ => vec![0; n_components],
    };
    let mut comps_open = comp_count
        .iter()
        .zip(&comp_target)
        .filter(|(n, t)| n > t)
        .count();
    let max_cost = match options.target {
        Target::MaxCost(c) => c,
        _ => f64::INFINITY,
    };

    let mut flatten_iterations = 0usize;
    let mut updates = Vec::with_capacity(32);
    let mut requeue: Vec<u32> = Vec::with_capacity(32);
    let mut t_removals = 0.0;
    let mut t_requeue = 0.0;
    let stop = loop {
        if run.mesh.n_vertices() <= total_target || comps_open == 0 {
            break StopReason::TargetReached;
        }
        let Some(Reverse(entry)) = run.heap.pop() else {
            break StopReason::QueueEmpty;
        };
        run.report.n_pops += 1;
        let v = entry.fine;
        if entry.version != run.version[v as usize] || !run.mesh.is_vertex(v) {
            run.report.n_stale_pops += 1;
            continue;
        }
        if entry.cost == f64::INFINITY {
            break StopReason::NothingRemovable;
        }
        if entry.cost > max_cost {
            break StopReason::MaxCost;
        }
        let comp = component[v as usize] as usize;
        if comp_count[comp] <= comp_target[comp] {
            run.report.n_component_done_pops += 1;
            continue;
        }

        // Commit (Alg. 1 lines 10-12).
        let t0 = Instant::now();
        let prepared: PreparedRemoval = match run.mesh.begin_removal(v) {
            Ok(p) => p,
            Err(f) => {
                *run.report
                    .commit_failures
                    .entry(refusal_name(&f.reason))
                    .or_default() += 1;
                run.last_refusal[v as usize] = Some(refusal_name(&f.reason));
                t_removals += t0.elapsed().as_secs_f64();
                continue;
            }
        };
        updates.clear();
        let cost = run.state.evaluate(
            &prepared,
            &options.metric,
            run.commit_frames,
            &mut run.alphas,
            Some(&mut updates),
        );
        let removal = run.mesh.finish_removal()?;
        run.n_removed += 1;
        // The transport (Eq. 5, 8).
        run.state.commit(&updates, v);
        let mut delaunay_flips: Vec<FlippedEdge> = Vec::new();
        if let Some(d) = &removal.delaunay {
            run.report.n_delaunay_refused += d.n_refused;
            delaunay_flips.extend_from_slice(&d.flips);
        }
        if compat.is_some_and(|r| r.delaunay_two_ring) {
            let seeds = two_ring_edges(run.mesh, removal.face as usize);
            let d = run.mesh.delaunay_pass(&seeds, true)?;
            run.report.n_delaunay_refused += d.n_refused;
            delaunay_flips.extend_from_slice(&d.flips);
        }
        t_removals += t0.elapsed().as_secs_f64();

        // Bookkeeping.
        match prepared.kind {
            RemovalKind::Interior => run.report.n_interior += 1,
            RemovalKind::Boundary => run.report.n_boundary += 1,
            RemovalKind::Ear => run.report.n_ear += 1,
        }
        comp_count[comp] -= 1;
        if comp_count[comp] == comp_target[comp] && comps_open > 0 {
            comps_open -= 1;
        }
        flatten_iterations += prepared.flatten.iterations;
        run.report.flatten_iterations_max = run
            .report
            .flatten_iterations_max
            .max(prepared.flatten.iterations);
        run.report.n_fallback_flips += prepared.fallback_flips.len();
        run.report.n_degree_flips +=
            prepared.degree_flips.len() + prepared.ear_flip.is_some() as usize;
        run.report.n_delaunay_flips += delaunay_flips.len();
        // S5: was the key computed after the last change of v and its neighbours?
        if options.check_keys {
            let ks = run.key_step[v as usize];
            let clean = run.touched[v as usize] <= ks
                && prepared
                    .neighbours
                    .iter()
                    .all(|n| run.touched[n.vertex as usize] <= ks);
            let k = &mut run.report.keys;
            if clean {
                k.n_clean += 1;
                if entry.cost.to_bits() != cost.to_bits() {
                    k.n_clean_mismatch += 1;
                }
            }
            if entry.cost.to_bits() != cost.to_bits() {
                k.n_differ += 1;
                let rel = (entry.cost - cost).abs() / cost.abs().max(f64::MIN_POSITIVE);
                k.max_relative_difference = k.max_relative_difference.max(rel);
            }
        }
        if options.record_sequence {
            run.report.sequence.push(RemovalRecord {
                vertex: v,
                kind: prepared.kind,
                key: entry.cost,
                cost,
            });
        }

        // The changed vertices (S6) and the re-queue.
        requeue.clear();
        requeue.extend(prepared.neighbours.iter().map(|n| n.vertex));
        let flips = prepared.flips().chain(delaunay_flips.iter());
        let mut flip_vertices: Vec<u32> = Vec::new();
        for f in flips {
            flip_vertices.extend_from_slice(&f.old);
            flip_vertices.extend_from_slice(&f.new);
        }
        run.touch(v);
        for &x in requeue.iter().chain(&flip_vertices) {
            run.touch(x);
        }
        if options.requeue == Requeue::StarsChanged {
            requeue.extend_from_slice(&flip_vertices);
            requeue.sort_unstable();
            requeue.dedup();
        }
        let t0 = Instant::now();
        for &x in &requeue {
            if x != v && run.mesh.is_vertex(x) {
                run.push(x);
            }
        }
        t_requeue += t0.elapsed().as_secs_f64();
        take_snapshots(&mut run, &mut previous_vertices);
    };

    let n_committed = run.report.n_interior + run.report.n_boundary + run.report.n_ear;
    let mut report = run.report;
    report.stop = stop;
    report.n_vertices_after = run.mesh.n_vertices();
    report.flatten_iterations_mean = if n_committed > 0 {
        flatten_iterations as f64 / n_committed as f64
    } else {
        0.0
    };
    report.n_clamped = run.mesh.n_clamped_points() - clamped_before;
    report.components = comp_before
        .iter()
        .zip(&comp_count)
        .map(|(&before, &after)| ComponentCount { before, after })
        .collect();
    for &v in run.mesh.vertices() {
        if let Some(name) = run.last_refusal[v as usize] {
            *report.skips.entry(name).or_default() += 1;
        }
    }
    report.timings.removals = t_removals;
    report.timings.requeue = t_requeue;
    report.channels = run.state;
    report.timings.total = t_start.elapsed().as_secs_f64();
    Ok(report)
}
