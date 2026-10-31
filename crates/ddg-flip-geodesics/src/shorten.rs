//! Geodesic paths, loops and networks, Alg. 2 `MakePathGeodesic` (SPEC §4.2–§4.4).
//!
//! # The queue
//!
//! A min-priority queue of joint sides keyed `(α, side, path, segment id)`, `α` the wedge
//! angle on that side and the segment the joint's outgoing one, compared with
//! `f64::total_cmp` and `Left < Right` (F4). A joint is queued with the sides
//! [`JointGeometry::sides_to_shorten`] gives (its smaller side if `α^min ≤ π − tol`, and the
//! larger one too if it is also below), never at a pinned vertex. Initially every joint of
//! every path, path by path in path order.
//!
//! # The loop (`gc: fg.cpp:919-975`)
//!
//! 1. Pop. **Stale** (F5) if the segment is gone, no longer has a predecessor, or the angle on
//!    that side recomputed now differs bitwise from the key: skip. Angles are summed corner by
//!    corner in a fixed order (F2), so an unchanged wedge gives the same bits.
//! 2. **Blocked** (F5): the wedge is not flexible: drop. It is queued again when a neighbour
//!    is replaced (below), or by the final sweep.
//! 3. If `stop_min_angle` is set and `α ≥` it: stop (`MinAngle`; the queue is ordered, so
//!    every live entry is at least as large). If `max_iterations` shortenings were made: stop
//!    (`IterationLimit`).
//! 4. Alg. 1 ([`crate::flipout`]), or the App. B move for a single-segment loop. One
//!    iteration. The running length is updated (`length += new − old`, F14).
//! 5. Re-queue (both sides, as above): the joint at the tail of every new segment, the joint
//!    at the tail of the segment that follows them (or that followed the removed run), and both
//!    joints of every segment that lay next to a removed segment in its edge's layer (the
//!    wedges at their ends may have become flexible; the reference re-queues the segment
//!    outermost on the wedge side, `gc: fg.cpp:175-182`). A joint that **stalled** (the outer
//!    arc was not shorter, F6) is remembered with its angle; it is not queued by a sweep again
//!    until that angle changes or a replacement re-queues it.
//! 6. If `stop_length_ratio` is set and the length fell below that fraction of the initial
//!    length: stop (`LengthRatio`).
//!
//! When the queue is empty, a **final sweep** recomputes every joint, in path order, and
//! queues every flexible side `≤ π − tol` that is not remembered as stalled (F5; the
//! reference stops at an empty queue). A sweep that queues nothing ends the run (`Geodesic`);
//! after `max_sweeps` sweeps that found work, the next one that finds work ends it
//! (`SweepLimit`). The queue is rebuilt without its stale entries when it holds more than
//! `4 m + 1024` entries (`m` segments); that never changes which entry is processed next.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BinaryHeap};
use std::f64::consts::PI;

use ddg_intrinsic::Triangulation;

use crate::error::FlipError;
use crate::flipout::{FlipOutResult, Outcome, flip_out};
use crate::network::{GeodesicNetwork, Joint, NONE, Network, SegmentId};
use crate::report::{ShortenOptions, ShortenReport, ShortenStatus, ShortenStep, Shortening};
use crate::scratch::{Located, Scratch, sides_to_shorten};
use crate::wedge::{JointGeometry, Side};

/// `log2` of the number of faces the corner-angle memo of a run holds
/// ([`crate::scratch::AngleCache`]) for a network of `m` segments: about `32 m`, between `2^9`
/// and `2^13`. Collisions cost time only.
fn angle_cache_bits(m: usize) -> u32 {
    (32 * m.max(1))
        .next_power_of_two()
        .trailing_zeros()
        .clamp(9, 13)
}

/// A queue entry: one side of the joint at the tail of `segment`.
#[derive(Debug, Clone, Copy)]
struct Entry {
    angle: f64,
    side: Side,
    path: u32,
    segment: u32,
}

impl Entry {
    fn key(&self) -> (Side, u32, u32) {
        (self.side, self.path, self.segment)
    }
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
        self.angle
            .total_cmp(&other.angle)
            .then_with(|| self.key().cmp(&other.key()))
    }
}

/// The joint at the tail of segment `id`, if it is live and has a predecessor.
fn joint_at(network: &Network, id: u32) -> Option<Joint> {
    let node = network.nodes.get(id as usize).copied().flatten()?;
    (node.prev != NONE).then_some(Joint {
        path: node.path as usize,
        incoming: SegmentId(node.prev),
        outgoing: SegmentId(id),
    })
}

/// Calls `f` with each side of the joint at the tail of `id` to queue
/// ([`JointGeometry::sides_to_shorten`]; none at a pinned vertex), with the joint located in
/// `scratch`.
fn for_each_side<T: Triangulation, F>(
    tri: &T,
    network: &Network,
    scratch: &mut Scratch,
    id: u32,
    tolerance: f64,
    mut f: F,
) -> Result<(), FlipError>
where
    F: FnMut(Entry, &Scratch, Located),
{
    let Some(joint) = joint_at(network, id) else {
        return Ok(());
    };
    if network.is_pinned(network.joint_vertex(tri.complex(), joint)) {
        return Ok(());
    }
    let at = scratch.locate(tri, network, joint)?;
    let (left, right) = scratch.angles(tri, at);
    let (sides, n) = sides_to_shorten(left, right, tolerance);
    for &(side, angle) in &sides[..n] {
        let entry = Entry {
            angle,
            side,
            path: joint.path as u32,
            segment: id,
        };
        f(entry, scratch, at);
    }
    Ok(())
}

/// Whether a popped entry is still current: its joint exists and the angle on its side is
/// bitwise the key. Returns the joint, located in `scratch`, if so.
fn current<T: Triangulation>(
    tri: &T,
    network: &Network,
    scratch: &mut Scratch,
    entry: &Entry,
) -> Result<Option<Located>, FlipError> {
    let Some(joint) = joint_at(network, entry.segment) else {
        return Ok(None);
    };
    let at = scratch.locate(tri, network, joint)?;
    let angle = scratch.angle(tri, at, entry.side);
    Ok((angle.to_bits() == entry.angle.to_bits()).then_some(at))
}

struct Run<'a> {
    options: &'a ShortenOptions,
    heap: BinaryHeap<Reverse<Entry>>,
    /// Joint sides that stalled, with the angle they had.
    stalled: BTreeMap<(u32, Side), u64>,
    report: ShortenReport,
    scratch: Scratch,
    /// The segments whose joints a replacement re-queues.
    ids: Vec<u32>,
}

/// A shortening whose queue is built ([`GeodesicNetwork::start_shorten`]) and whose loop has
/// not run yet.
pub struct ShortenSession<'n, T: Triangulation> {
    network: &'n mut GeodesicNetwork<T>,
    run: Run<'n>,
    max_iterations: usize,
}

impl<T: Triangulation> ShortenSession<'_, T> {
    /// Runs the rest of [`GeodesicNetwork::shorten`].
    ///
    /// # Errors
    /// As [`GeodesicNetwork::shorten`].
    pub fn run(self) -> Result<ShortenReport, FlipError> {
        self.finish::<fn(&ShortenStep, &GeodesicNetwork<T>)>(None)
    }

    /// Runs the rest of [`GeodesicNetwork::shorten_with`].
    ///
    /// # Errors
    /// As [`GeodesicNetwork::shorten`].
    pub fn run_with<F>(self, mut observer: F) -> Result<ShortenReport, FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        self.finish(Some(&mut observer))
    }

    fn finish<F>(self, observer: Option<&mut F>) -> Result<ShortenReport, FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        let ShortenSession {
            network,
            run,
            max_iterations,
        } = self;
        network.shorten_loop(run, max_iterations, observer)
    }
}

impl<T: Triangulation> GeodesicNetwork<T> {
    /// Makes every path of the network geodesic: Alg. 2 of the paper with the decisions of
    /// SPEC §4.2–§4.4 (see the [module documentation](crate::shorten)). Paths, loops and
    /// networks alike; joints at pinned vertices (and the ends of open paths) are never
    /// shortened.
    ///
    /// The result is not certified here: [`certify`](Self::certify) does that on the final
    /// state (SPEC §7 layer 1a).
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] (nothing is changed then); [`FlipError::Internal`] for a
    /// broken invariant (a bug); with [`ShortenOptions::check_invariants`], the errors of
    /// [`Network::check`] after the iteration that broke it.
    pub fn shorten(&mut self, options: &ShortenOptions) -> Result<ShortenReport, FlipError> {
        self.run_shorten::<fn(&ShortenStep, &GeodesicNetwork<T>)>(options, None)
    }

    /// [`shorten`](Self::shorten), calling `observer` after every iteration with what it did
    /// and the network as it is then (for tests: monotone length, isotopy, the Δ-complex
    /// cases).
    ///
    /// # Errors
    /// As [`shorten`](Self::shorten).
    pub fn shorten_with<F>(
        &mut self,
        options: &ShortenOptions,
        mut observer: F,
    ) -> Result<ShortenReport, FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        self.run_shorten(options, Some(&mut observer))
    }

    /// The shortening, with or without an observer (without one, no [`ShortenStep`] is
    /// built).
    fn run_shorten<F>(
        &mut self,
        options: &ShortenOptions,
        observer: Option<&mut F>,
    ) -> Result<ShortenReport, FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        self.start_shorten(options)?.finish(observer)
    }

    /// The first part of [`shorten`](Self::shorten): checks the options and queues every joint
    /// of every path; [`ShortenSession::run`] does the rest. `start_shorten(o)?.run()` is
    /// `shorten(o)`, result for result. For timing the loop alone, as the reference does (it
    /// builds its queue when the network is made, `gc: fg.cpp:47`).
    ///
    /// # Errors
    /// As [`shorten`](Self::shorten) (nothing is changed).
    pub fn start_shorten<'n>(
        &'n mut self,
        options: &'n ShortenOptions,
    ) -> Result<ShortenSession<'n, T>, FlipError> {
        options.validate()?;
        let tol = options.straight_tolerance;
        let m0 = self.network.n_segments();
        let max_iterations = options
            .max_iterations
            .unwrap_or(100 * (self.tri.complex().n_edges() + m0));
        let initial_length = self.length;
        let mut run = Run {
            options,
            heap: BinaryHeap::new(),
            stalled: BTreeMap::new(),
            report: ShortenReport {
                status: ShortenStatus::Geodesic,
                iterations: 0,
                flips: 0,
                refused_flips: 0,
                stalled: 0,
                sweeps: 0,
                loop_perturbations: 0,
                contractions: 0,
                twice_incident_wedges: 0,
                degree_one_wedges: 0,
                stale_entries: 0,
                blocked_entries: 0,
                contracted: Vec::new(),
                initial_length,
                final_length: initial_length,
            },
            scratch: Scratch::new(angle_cache_bits(m0)),
            ids: Vec::new(),
        };
        for p in 0..self.network.n_paths() {
            for joint in self.network.path_joints(p) {
                let heap = &mut run.heap;
                for_each_side(
                    &self.tri,
                    &self.network,
                    &mut run.scratch,
                    joint.outgoing.0,
                    tol,
                    |entry, _, _| heap.push(Reverse(entry)),
                )?;
            }
        }
        Ok(ShortenSession {
            network: self,
            run,
            max_iterations,
        })
    }

    /// The loop of Alg. 2 on a queue built by [`start_shorten`](Self::start_shorten).
    fn shorten_loop<F>(
        &mut self,
        mut run: Run<'_>,
        max_iterations: usize,
        mut observer: Option<&mut F>,
    ) -> Result<ShortenReport, FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        let options = run.options;
        let tol = options.straight_tolerance;
        let initial_length = run.report.initial_length;
        let status = 'run: loop {
            while let Some(Reverse(entry)) = run.heap.pop() {
                let Some(at) = current(&self.tri, &self.network, &mut run.scratch, &entry)? else {
                    run.report.stale_entries += 1;
                    continue;
                };
                if self.network.is_pinned(at.vertex) {
                    continue;
                }
                if !run
                    .scratch
                    .is_flexible(&self.network, &self.tri, at, entry.side)
                {
                    run.report.blocked_entries += 1;
                    continue;
                }
                if run.stalled.get(&(entry.segment, entry.side)) == Some(&entry.angle.to_bits()) {
                    // A second entry of a joint side that stalled with this very angle.
                    continue;
                }
                if let Some(m) = options.stop_min_angle
                    && entry.angle >= m
                {
                    break 'run ShortenStatus::MinAngle;
                }
                if run.report.iterations >= max_iterations {
                    break 'run ShortenStatus::IterationLimit;
                }
                self.step(&mut run, at, entry, &mut observer)?;
                if let Some(r) = options.stop_length_ratio
                    && self.length < r * initial_length
                {
                    break 'run ShortenStatus::LengthRatio;
                }
                if run.heap.len() > 4 * self.network.n_segments() + 1024 {
                    let entries = std::mem::take(&mut run.heap).into_vec();
                    for Reverse(e) in entries {
                        if current(&self.tri, &self.network, &mut run.scratch, &e)?.is_some() {
                            run.heap.push(Reverse(e));
                        } else {
                            run.report.stale_entries += 1;
                        }
                    }
                }
            }
            // The final sweep.
            let mut found = Vec::new();
            for p in 0..self.network.n_paths() {
                for joint in self.network.path_joints(p) {
                    let stalled = &run.stalled;
                    let (network, tri) = (&self.network, &self.tri);
                    for_each_side(
                        &self.tri,
                        &self.network,
                        &mut run.scratch,
                        joint.outgoing.0,
                        tol,
                        |entry, scratch, at| {
                            let is_stalled = stalled.get(&(entry.segment, entry.side))
                                == Some(&entry.angle.to_bits());
                            if scratch.is_flexible(network, tri, at, entry.side) && !is_stalled {
                                found.push(entry);
                            }
                        },
                    )?;
                }
            }
            if found.is_empty() {
                break 'run ShortenStatus::Geodesic;
            }
            if run.report.sweeps >= options.max_sweeps {
                break 'run ShortenStatus::SweepLimit;
            }
            run.report.sweeps += 1;
            run.heap.extend(found.into_iter().map(Reverse));
        };
        run.report.status = status;
        run.report.final_length = self.length;
        run.report.contracted.sort_unstable();
        Ok(run.report)
    }

    /// One iteration: Alg. 1 on the popped joint side (`at`, located last in the scratch), the
    /// bookkeeping, the observer.
    fn step<F>(
        &mut self,
        run: &mut Run<'_>,
        at: Located,
        entry: Entry,
        observer: &mut Option<&mut F>,
    ) -> Result<(), FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        let joint = at.joint;
        let path = joint.path;
        let result = flip_out(
            &mut self.tri,
            &mut self.network,
            at,
            entry.side,
            &mut run.scratch,
        )?;
        let tol = run.options.straight_tolerance;
        let report = &mut run.report;
        report.iterations += 1;
        report.flips += result.stats.flips;
        report.refused_flips += result.stats.refused;
        report.twice_incident_wedges += usize::from(result.stats.twice_incident);
        report.degree_one_wedges += usize::from(result.stats.degree_one);
        match result.outcome {
            Outcome::Stalled { .. } => {
                report.stalled += 1;
                // Remember the joint side with its angle after the flips.
                let now = run.scratch.locate(&self.tri, &self.network, joint)?;
                let angle = run.scratch.angle(&self.tri, now, entry.side);
                run.stalled
                    .insert((entry.segment, entry.side), angle.to_bits());
            }
            Outcome::LoopMoved { .. } => report.loop_perturbations += 1,
            Outcome::Contracted {
                path_contracted, ..
            } => {
                report.contractions += 1;
                if path_contracted {
                    report.contracted.push(path);
                }
            }
            Outcome::Replaced { .. } => {}
        }
        self.update_length(result.outcome);
        if result.outcome.replaced() {
            let replaced = &run.scratch.replaced;
            run.ids.clear();
            run.ids.extend_from_slice(&replaced.new_ids);
            run.ids.extend(replaced.after);
            for &n in &replaced.layer_neighbours {
                run.ids.push(n);
                if let Some(node) = self.network.nodes.get(n as usize).copied().flatten()
                    && node.next != NONE
                {
                    run.ids.push(node.next);
                }
            }
            for &id in &run.ids {
                run.stalled.remove(&(id, Side::Left));
                run.stalled.remove(&(id, Side::Right));
                let heap = &mut run.heap;
                for_each_side(
                    &self.tri,
                    &self.network,
                    &mut run.scratch,
                    id,
                    tol,
                    |e, _, _| heap.push(Reverse(e)),
                )?;
            }
        }
        if run.options.check_invariants {
            self.check_invariants()?;
        }
        if let Some(observer) = observer {
            let FlipOutResult { outcome, stats } = result;
            let step = ShortenStep {
                iteration: run.report.iterations,
                joint,
                side: entry.side,
                vertex: at.vertex,
                angle: entry.angle,
                flips: stats.flips,
                refused_flips: stats.refused,
                shortening: outcome.shortening(&run.scratch.replaced.new_ids),
                length: self.length,
            };
            observer(&step, self);
        }
        Ok(())
    }

    /// The running length (F14): `length += new − old` per step; when a path contracts to
    /// nothing the length is summed afresh, so that a network of contracted paths has length
    /// exactly `0` rather than a rounding residue.
    fn update_length(&mut self, outcome: Outcome) {
        match outcome {
            Outcome::Contracted {
                path_contracted: true,
                ..
            } => self.length = self.recompute_length(),
            _ => self.length += outcome.length_change(),
        }
    }

    /// Shortens one joint on one side by Alg. 1 (or the App. B move of a single-segment
    /// loop), as one iteration of [`shorten`](Self::shorten) does, and updates the running
    /// length. Nothing else is queued or swept.
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] if `joint` is not two consecutive live segments;
    /// [`FlipError::NotShortenable`] (nothing changed) if its vertex is pinned, the wedge on
    /// `side` crosses the boundary, is not flexible, or its angle is not below `π`;
    /// [`FlipError::Internal`] for a broken invariant.
    pub fn shorten_joint(&mut self, joint: Joint, side: Side) -> Result<Shortening, FlipError> {
        let geometry = JointGeometry::new(&self.tri, &self.network, joint)?;
        let vertex = geometry.vertex();
        let refuse = |reason| {
            Err(FlipError::NotShortenable {
                vertex,
                side,
                reason,
            })
        };
        if self.network.is_pinned(vertex) {
            return refuse("the vertex is pinned");
        }
        let wedge = geometry.wedge(side);
        if wedge.boundary {
            return refuse("the wedge crosses the boundary");
        }
        if !geometry.is_flexible(side) {
            return refuse("the wedge is not flexible (it contains a segment)");
        }
        if wedge.angle.partial_cmp(&PI) != Some(std::cmp::Ordering::Less) {
            return refuse("the wedge angle is not below pi");
        }
        let mut scratch = Scratch::new(6);
        let at = scratch.locate(&self.tri, &self.network, joint)?;
        let result = flip_out(&mut self.tri, &mut self.network, at, side, &mut scratch)?;
        self.update_length(result.outcome);
        Ok(result.outcome.shortening(&scratch.replaced.new_ids))
    }
}
