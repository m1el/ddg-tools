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
use crate::flipout::{FlipOutResult, flip_out};
use crate::network::{GeodesicNetwork, Joint, NONE, Network, SegmentId};
use crate::report::{ShortenOptions, ShortenReport, ShortenStatus, ShortenStep, Shortening};
use crate::wedge::{JointGeometry, Side};

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

/// The sides of the joint at the tail of `id` to queue (none at a pinned vertex).
fn entries_for<T: Triangulation>(
    tri: &T,
    network: &Network,
    id: u32,
    tolerance: f64,
) -> Result<Vec<(Entry, JointGeometry)>, FlipError> {
    let Some(joint) = joint_at(network, id) else {
        return Ok(Vec::new());
    };
    if network.is_pinned(network.joint_vertex(tri.complex(), joint)) {
        return Ok(Vec::new());
    }
    let geometry = JointGeometry::new(tri, network, joint)?;
    Ok(geometry
        .sides_to_shorten(tolerance)
        .into_iter()
        .map(|(side, angle)| {
            (
                Entry {
                    angle,
                    side,
                    path: joint.path as u32,
                    segment: id,
                },
                geometry.clone(),
            )
        })
        .collect())
}

/// Whether a popped entry is still current: its joint exists and the angle on its side is
/// bitwise the key. Returns the joint's geometry if so.
fn current<T: Triangulation>(
    tri: &T,
    network: &Network,
    entry: &Entry,
) -> Result<Option<JointGeometry>, FlipError> {
    let Some(joint) = joint_at(network, entry.segment) else {
        return Ok(None);
    };
    let geometry = JointGeometry::new(tri, network, joint)?;
    Ok((geometry.angle(entry.side).to_bits() == entry.angle.to_bits()).then_some(geometry))
}

struct Run<'a> {
    options: &'a ShortenOptions,
    heap: BinaryHeap<Reverse<Entry>>,
    /// Joint sides that stalled, with the angle they had.
    stalled: BTreeMap<(u32, Side), u64>,
    report: ShortenReport,
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
        self.shorten_with(options, |_, _| {})
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
        };
        for p in 0..self.network.n_paths() {
            for joint in self.network.path_joints(p) {
                for (entry, _) in entries_for(&self.tri, &self.network, joint.outgoing.0, tol)? {
                    run.heap.push(Reverse(entry));
                }
            }
        }

        let status = 'run: loop {
            while let Some(Reverse(entry)) = run.heap.pop() {
                let Some(geometry) = current(&self.tri, &self.network, &entry)? else {
                    run.report.stale_entries += 1;
                    continue;
                };
                if self.network.is_pinned(geometry.vertex()) {
                    continue;
                }
                if !geometry.is_flexible(entry.side) {
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
                self.step(&mut run, &geometry, entry, &mut observer)?;
                if let Some(r) = options.stop_length_ratio
                    && self.length < r * initial_length
                {
                    break 'run ShortenStatus::LengthRatio;
                }
                if run.heap.len() > 4 * self.network.n_segments() + 1024 {
                    let entries = std::mem::take(&mut run.heap).into_vec();
                    for Reverse(e) in entries {
                        if current(&self.tri, &self.network, &e)?.is_some() {
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
                    for (entry, geometry) in
                        entries_for(&self.tri, &self.network, joint.outgoing.0, tol)?
                    {
                        let stalled = run.stalled.get(&(entry.segment, entry.side))
                            == Some(&entry.angle.to_bits());
                        if geometry.is_flexible(entry.side) && !stalled {
                            found.push(entry);
                        }
                    }
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

    /// One iteration: Alg. 1 on the popped joint side, the bookkeeping, the observer.
    fn step<F>(
        &mut self,
        run: &mut Run<'_>,
        geometry: &JointGeometry,
        entry: Entry,
        observer: &mut F,
    ) -> Result<(), FlipError>
    where
        F: FnMut(&ShortenStep, &GeodesicNetwork<T>),
    {
        let joint = geometry.joint;
        let path = joint.path;
        let result = flip_out(&mut self.tri, &mut self.network, joint, entry.side)?;
        let tol = run.options.straight_tolerance;
        let report = &mut run.report;
        report.iterations += 1;
        report.flips += result.stats.flips;
        report.refused_flips += result.stats.refused;
        report.twice_incident_wedges += usize::from(result.stats.twice_incident);
        report.degree_one_wedges += usize::from(result.stats.degree_one);
        match &result.shortening {
            Shortening::Stalled { .. } => {
                report.stalled += 1;
                // Remember the joint side with its angle after the flips.
                let now = JointGeometry::new(&self.tri, &self.network, joint)?;
                run.stalled
                    .insert((entry.segment, entry.side), now.angle(entry.side).to_bits());
            }
            Shortening::LoopMoved { .. } => report.loop_perturbations += 1,
            Shortening::Contracted {
                path_contracted, ..
            } => {
                report.contractions += 1;
                if *path_contracted {
                    report.contracted.push(path);
                }
            }
            Shortening::Replaced { .. } => {}
        }
        self.update_length(&result.shortening);
        if let Some(replaced) = &result.replaced {
            let mut ids: Vec<u32> = replaced.new_ids.clone();
            ids.extend(replaced.after);
            for &n in &replaced.layer_neighbours {
                ids.push(n);
                if let Some(node) = self.network.nodes.get(n as usize).copied().flatten()
                    && node.next != NONE
                {
                    ids.push(node.next);
                }
            }
            for id in ids {
                run.stalled.remove(&(id, Side::Left));
                run.stalled.remove(&(id, Side::Right));
                for (e, _) in entries_for(&self.tri, &self.network, id, tol)? {
                    run.heap.push(Reverse(e));
                }
            }
        }
        if run.options.check_invariants {
            self.check_invariants()?;
        }
        let FlipOutResult {
            shortening, stats, ..
        } = result;
        let step = ShortenStep {
            iteration: run.report.iterations,
            joint,
            side: entry.side,
            vertex: geometry.vertex(),
            angle: entry.angle,
            flips: stats.flips,
            refused_flips: stats.refused,
            shortening,
            length: self.length,
        };
        observer(&step, self);
        Ok(())
    }

    /// The running length (F14): `length += new − old` per step; when a path contracts to
    /// nothing the length is summed afresh, so that a network of contracted paths has length
    /// exactly `0` rather than a rounding residue.
    fn update_length(&mut self, shortening: &Shortening) {
        match shortening {
            Shortening::Contracted {
                path_contracted: true,
                ..
            } => self.length = self.recompute_length(),
            _ => self.length += shortening.length_change(),
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
        let result = flip_out(&mut self.tri, &mut self.network, joint, side)?;
        self.update_length(&result.shortening);
        Ok(result.shortening)
    }
}
