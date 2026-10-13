//! Intrinsic Delaunay refinement (paper Sec. 4.2): vertices are inserted until no triangle
//! has an angle below a bound, while the triangulation stays Delaunay.
//!
//! The algorithm is the planar one carried to the surface. On a closed surface without sharp
//! cone points it is the paper's: flip to Delaunay, take a triangle with an angle below the
//! bound, insert its circumcentre, repeat. The circumcentre is found by tracing from the
//! triangle over the triangulation. Two things the paper leaves open are decided here and
//! written out in SPEC N1 and N4:
//!
//! **Boundary (N1).** Ruppert's rule. A boundary edge is *encroached* when the angle opposite
//! it exceeds 90° (a vertex lies inside the circle on the edge as diameter). Encroached
//! boundary edges are split first. A circumcentre is not inserted when the path to it leaves
//! the surface, or when it would encroach a boundary edge of the triangles its insertion
//! would remove; that boundary edge is split instead. A boundary edge is split at its
//! midpoint, except next to a boundary vertex with an interior angle below 60°, where it is
//! split at a distance from that vertex that is a power of two times a fixed length (so that
//! the two boundary edges at the vertex are cut at the same distances and stop encroaching on
//! each other).
//!
//! **Small cone angles (N4).** A corner at a vertex whose angle sum is below the bound can
//! never meet the bound; it is exempt. A circumcentre next to an interior vertex with an
//! angle sum below 60° would get an edge to itself around that vertex that is shorter than
//! every edge of the triangle it was meant to improve; it is not inserted (this contains the
//! paper's rule, which is the case of a vertex of degree one).
//!
//! Everything that is left below the bound is listed in the report with the reason.

use std::collections::{BinaryHeap, HashMap};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, TAU};

use ddg_core::topology::NONE;

use crate::delaunay::{DelaunayCriterion, EdgeQueue};
use crate::error::{InsertRefusal, IntrinsicError};
use crate::formulas::{circumcenter_barycentric, triangle_layout, unfold_apex};
use crate::insert::Insertion;
use crate::point::SurfacePoint;
use crate::trace::{TraceOptions, TraceStop};
use crate::triangulation::IntrinsicTriangulation;

/// How many vertices [`IntrinsicTriangulation::refine`] may insert.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsertionLimit {
    /// At most this many for every vertex the triangulation has when refinement starts, and
    /// at least 1000 in all (a small mesh with thin triangles can need many times its own
    /// number of vertices).
    PerVertex(usize),
    /// At most this many.
    Vertices(usize),
    /// No limit.
    Unlimited,
}

/// Options of [`IntrinsicTriangulation::refine`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineOptions {
    /// The angle bound in radians: refinement ends when no triangle has a corner below it,
    /// except those listed in the report. Between 0 and 35°; with 0 only encroached boundary
    /// edges are split.
    pub min_angle: f64,
    /// If given, a triangle with a larger circumradius is refined as well.
    pub max_circumradius: Option<f64>,
    /// The limit on inserted vertices; reaching it is reported, not an error.
    pub max_insertions: InsertionLimit,
    /// The Delaunay test used for the flips.
    pub criterion: DelaunayCriterion,
}

impl Default for RefineOptions {
    /// A bound of 25°, no bound on the circumradius, at most 10 insertions per vertex, the
    /// default Delaunay test. 25° is the bound up to which the planar algorithm is known to
    /// end (26.5° for Chew's second algorithm); the paper's 30° is reached in practice too and
    /// costs about two and a half times as many vertices (the crate's README has the table).
    fn default() -> Self {
        RefineOptions {
            min_angle: 25f64.to_radians(),
            max_circumradius: None,
            max_insertions: InsertionLimit::PerVertex(10),
            criterion: DelaunayCriterion::default(),
        }
    }
}

/// How refinement ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefineStatus {
    /// No triangle is left to refine: every corner meets the bound except those listed in
    /// [`RefineReport::below_bound`], each with its reason.
    Done,
    /// The limit on inserted vertices was reached.
    InsertionLimit,
}

/// Why a corner is below the angle bound after refinement.
#[derive(Debug, Clone, PartialEq)]
pub enum BelowReason {
    /// The corner is at a vertex whose angle sum is below the bound: no triangulation has a
    /// corner of the bound there.
    AngleSum {
        /// The angle sum of the vertex.
        angle_sum: f64,
    },
    /// The circumcentre of the triangle would lie so close to an interior vertex with an
    /// angle sum below 60° that the new vertex would get an edge to itself around that vertex
    /// shorter than the shortest edge of the triangle.
    NeedleVertex {
        /// The vertex with the small angle sum.
        vertex: u32,
        /// Its angle sum.
        angle_sum: f64,
    },
    /// The circumcentre could not be inserted (the insertion or the edge split was refused,
    /// or the path to it ran exactly into a vertex).
    NotInsertable {
        /// What was refused.
        detail: String,
    },
    /// The limit on inserted vertices was reached before the triangle's turn.
    InsertionLimit,
}

/// A corner that is below the angle bound after refinement.
#[derive(Debug, Clone, PartialEq)]
pub struct BelowBound {
    /// The face.
    pub face: u32,
    /// The vertex of the corner.
    pub vertex: u32,
    /// The angle of the corner.
    pub angle: f64,
    /// Why it was left.
    pub reason: BelowReason,
}

/// What [`IntrinsicTriangulation::refine`] did and what it left.
#[derive(Debug, Clone, PartialEq)]
pub struct RefineReport {
    /// How it ended.
    pub status: RefineStatus,
    /// Number of vertices inserted: the sum of the next three numbers.
    pub n_inserted: usize,
    /// Circumcentres inserted inside a face.
    pub n_face_insertions: usize,
    /// Interior edges split because a circumcentre lay on them (within the shape limit of an
    /// insertion, SPEC N5).
    pub n_edge_splits: usize,
    /// Boundary edges split.
    pub n_boundary_splits: usize,
    /// Number of flips, including those of the first flip to Delaunay.
    pub n_flips: usize,
    /// Number of new vertices whose traced place was moved onto an input edge
    /// ([`LOCATION_SNAP`](crate::insert::LOCATION_SNAP)).
    pub n_snapped: usize,
    /// The smallest corner angle of the result.
    pub min_angle: f64,
    /// The largest corner angle of the result.
    pub max_angle: f64,
    /// The smallest corner angle over the corners at vertices whose angle sum is not below
    /// the bound (`+∞` if there is none): the smallest angle that a triangulation could
    /// improve. It is at least the bound unless [`below_bound`](RefineReport::below_bound)
    /// lists a corner for another reason than its vertex's angle sum.
    pub min_angle_reachable: f64,
    /// Every corner below the bound, with the reason, in the order of the faces.
    pub below_bound: Vec<BelowBound>,
    /// Faces whose circumradius exceeds [`RefineOptions::max_circumradius`] at the end.
    pub too_large: Vec<u32>,
    /// Interior edges that fail the Delaunay test at the end (zero unless a flip was refused).
    pub n_not_delaunay: usize,
    /// Boundary edges that are encroached at the end (zero unless a split was refused or the
    /// limit was reached).
    pub n_encroached: usize,
}

/// What the triangles around a candidate point say about inserting it.
enum Verdict {
    Clear,
    Encroaches(u32),
    Needle(u32),
}

struct Refiner<'a> {
    tri: &'a mut IntrinsicTriangulation,
    options: RefineOptions,
    queue: EdgeQueue,
    /// Faces to look at: (area, face, stamp of the face when it was put here).
    heap: BinaryHeap<(u64, u32, u32)>,
    stamps: Vec<u32>,
    /// Boundary edges to split.
    encroached: Vec<u32>,
    /// Faces that cannot be refined, with the stamp they had, and why.
    blocked: HashMap<u32, (u32, BelowReason)>,
    /// Boundary vertices with an interior angle below 60°: the length whose multiples by
    /// powers of two are the distances at which their boundary edges are split.
    shells: HashMap<u32, f64>,
    touched: Vec<u32>,
    /// Per face: the number of the search that visited it last (see `examine`).
    mark: Vec<u32>,
    search: u32,
    /// Per face: how often it was put back unchanged after a boundary split made for it.
    retries: HashMap<u32, (u32, u32)>,
    budget: usize,
    n_face_insertions: usize,
    n_edge_splits: usize,
    n_boundary_splits: usize,
    n_flips: usize,
    n_snapped: usize,
}

/// How far an angle opposite a boundary edge may exceed 90° before the edge counts as
/// encroached.
const ENCROACHMENT_TOLERANCE: f64 = 1e-9;

impl Refiner<'_> {
    fn n_inserted(&self) -> usize {
        self.n_face_insertions + self.n_edge_splits + self.n_boundary_splits
    }

    /// Whether a corner at vertex `v` can never meet the bound.
    fn exempt(&self, v: usize) -> bool {
        self.tri.angle_sums[v] < self.options.min_angle
    }

    /// The smallest angle of face `f` among the corners that are not exempt, with its
    /// vertex, if it is below the bound.
    fn small_corner(&self, f: usize) -> Option<(f64, u32)> {
        let mut worst: Option<(f64, u32)> = None;
        for k in 0..3 {
            let v = self.tri.faces[f][k];
            if self.exempt(v as usize) {
                continue;
            }
            let angle = self.tri.corner_angle(3 * f + k);
            if angle < self.options.min_angle && worst.is_none_or(|(a, _)| angle < a) {
                worst = Some((angle, v));
            }
        }
        worst
    }

    fn circumradius(&self, f: usize) -> f64 {
        let [a, b, c] = self.tri.face_lengths(f);
        a * b * c / (4.0 * self.tri.face_area(f))
    }

    fn too_large(&self, f: usize) -> bool {
        self.options
            .max_circumradius
            .is_some_and(|r| self.circumradius(f) > r)
    }

    fn needs_refinement(&self, f: usize) -> bool {
        self.small_corner(f).is_some() || self.too_large(f)
    }

    fn is_encroached(&self, e: usize) -> bool {
        let h = self.tri.edge_halfedge[e] as usize;
        self.tri.twin[h] == NONE
            && self.tri.corner_angle(self.tri.prev(h)) > FRAC_PI_2 + ENCROACHMENT_TOLERANCE
    }

    /// Notes that face `f` changed: it is looked at again.
    fn touch(&mut self, f: usize) {
        if f >= self.stamps.len() {
            self.stamps.resize(f + 1, 0);
        }
        self.stamps[f] += 1;
        if self.needs_refinement(f) {
            self.heap
                .push((self.tri.face_area(f).to_bits(), f as u32, self.stamps[f]));
        }
        for k in 0..3 {
            let h = 3 * f + k;
            if self.tri.twin[h] == NONE {
                let e = self.tri.halfedge_edge[h];
                if self.is_encroached(e as usize) {
                    self.encroached.push(e);
                }
            }
        }
    }

    /// Flips back to Delaunay from the queue and looks at every face that changed.
    fn restore_delaunay(&mut self) {
        self.touched.clear();
        let run = self.tri.run_flip_queue(
            &mut self.queue,
            &self.options.criterion,
            usize::MAX,
            Some(&mut self.touched),
        );
        self.n_flips += run.n_flips;
        let mut touched = std::mem::take(&mut self.touched);
        touched.sort_unstable();
        touched.dedup();
        for &f in &touched {
            self.touch(f as usize);
        }
        self.touched = touched;
    }

    /// After a vertex was inserted: its faces are new, their edges may have to be flipped.
    fn after_insertion(&mut self, insertion: &Insertion) {
        self.n_snapped += usize::from(insertion.snapped.is_some());
        let v = insertion.vertex as usize;
        let faces: Vec<usize> = self.tri.outgoing_halfedges(v).map(|h| h / 3).collect();
        for &f in &faces {
            for k in 0..3 {
                self.queue.push(self.tri.halfedge_edge[3 * f + k]);
            }
        }
        for &f in &faces {
            self.touch(f);
        }
        self.restore_delaunay();
    }

    /// The parameter at which boundary edge `e` is split: its midpoint, or next to a boundary
    /// vertex with a small angle the nearest power-of-two distance.
    fn boundary_split_parameter(&self, e: usize) -> f64 {
        let [a, b] = self.tri.edge_vertices(e);
        let length = self.tri.edge_lengths[e];
        let shell = |unit: f64| {
            let exponent = (0.5 * length / unit).log2().round();
            (unit * exponent.exp2() / length).clamp(0.25, 0.75)
        };
        match (self.shells.get(&(a as u32)), self.shells.get(&(b as u32))) {
            (Some(&unit), None) => shell(unit),
            (None, Some(&unit)) => 1.0 - shell(unit),
            _ => 0.5,
        }
    }

    /// Splits boundary edge `e`. Returns whether it was split.
    fn split_boundary(&mut self, e: usize) -> Result<bool, IntrinsicError> {
        let t = self.boundary_split_parameter(e);
        match self.tri.split_edge(e, t) {
            Ok(insertion) => {
                self.n_boundary_splits += 1;
                self.after_insertion(&insertion);
                Ok(true)
            }
            Err(IntrinsicError::InsertionRefused { .. }) => Ok(false),
            Err(other) => Err(other),
        }
    }

    /// Looks at the triangles that the insertion of the point `bary` of face `g` would
    /// remove: those whose circumcircle contains it, unfolded around it.
    fn examine(&mut self, g: usize, bary: [f64; 3], shortest: f64) -> Verdict {
        self.search += 1;
        let search = self.search;
        self.mark.resize(self.tri.n_faces(), 0);
        let mark = &mut self.mark;
        let tri = &*self.tri;
        let Some(v) = triangle_layout(tri.face_lengths(g)) else {
            return Verdict::Clear;
        };
        let c = [
            bary[1] * v[1][0] + bary[2] * v[2][0],
            bary[1] * v[1][1] + bary[2] * v[2][1],
        ];
        mark[g] = search;
        let mut stack: Vec<(usize, [[f64; 2]; 3])> = vec![(g, v)];
        let mut encroached: Option<u32> = None;
        while let Some((f, q)) = stack.pop() {
            for k in 0..3 {
                let h = 3 * f + k;
                let (a, b) = (q[k], q[(k + 1) % 3]);
                // An interior vertex with a small angle sum: the loop around it.
                let vertex = tri.faces[f][k] as usize;
                let theta = tri.angle_sums[vertex];
                if theta < FRAC_PI_3 && !tri.is_boundary_vertex(vertex) {
                    let d = (a[0] - c[0]).hypot(a[1] - c[1]);
                    if 2.0 * d * (0.5 * theta).sin() < shortest * (1.0 - 1e-9) {
                        return Verdict::Needle(vertex as u32);
                    }
                }
                let t = tri.twin[h];
                if t == NONE {
                    // Inside the circle on the boundary edge as diameter: the angle at the
                    // point is obtuse.
                    let dot = (a[0] - c[0]) * (b[0] - c[0]) + (a[1] - c[1]) * (b[1] - c[1]);
                    if dot < 0.0 && encroached.is_none() {
                        encroached = Some(tri.halfedge_edge[h]);
                    }
                    continue;
                }
                let t = t as usize;
                let (nf, kt) = (t / 3, t % 3);
                if mark[nf] == search {
                    continue;
                }
                let l = tri.face_lengths(nf);
                let Some(apex) = unfold_apex(b, a, l[kt], l[(kt + 1) % 3], l[(kt + 2) % 3]) else {
                    continue;
                };
                // Whether the point is inside the circumcircle of (b, a, apex).
                let rel = [b, a, apex].map(|p| {
                    let d = [p[0] - c[0], p[1] - c[1]];
                    [d[0], d[1], d[0] * d[0] + d[1] * d[1]]
                });
                let det = rel[0][0] * (rel[1][1] * rel[2][2] - rel[1][2] * rel[2][1])
                    - rel[0][1] * (rel[1][0] * rel[2][2] - rel[1][2] * rel[2][0])
                    + rel[0][2] * (rel[1][0] * rel[2][1] - rel[1][1] * rel[2][0]);
                if det > 0.0 {
                    mark[nf] = search;
                    let mut q2 = [[0.0; 2]; 3];
                    q2[kt] = b;
                    q2[(kt + 1) % 3] = a;
                    q2[(kt + 2) % 3] = apex;
                    stack.push((nf, q2));
                }
            }
        }
        match encroached {
            Some(e) => Verdict::Encroaches(e),
            None => Verdict::Clear,
        }
    }

    fn block(&mut self, f: usize, reason: BelowReason) {
        self.blocked.insert(f as u32, (self.stamps[f], reason));
    }

    /// One step for face `f`, which needs refinement.
    fn refine_face(&mut self, f: usize) -> Result<(), IntrinsicError> {
        let l = self.tri.face_lengths(f);
        let shortest = l[0].min(l[1]).min(l[2]);
        let (Some(center), Some(layout)) = (circumcenter_barycentric(l), triangle_layout(l)) else {
            self.block(
                f,
                BelowReason::NotInsertable {
                    detail: "the circumcentre of the face cannot be computed".to_string(),
                },
            );
            return Ok(());
        };
        // Where the circumcentre is: in the face itself, or at the end of the path from its
        // barycentre.
        let (g, bary) = if center.iter().all(|&b| b > 0.0) {
            (f, center)
        } else {
            let third = 1.0 / 3.0;
            let p = [0, 1].map(|i| {
                center[0] * layout[0][i] + center[1] * layout[1][i] + center[2] * layout[2][i]
            });
            let o = [0, 1].map(|i| third * (layout[0][i] + layout[1][i] + layout[2][i]));
            let d = [p[0] - o[0], p[1] - o[1]];
            let direction = d[1].atan2(d[0]).rem_euclid(TAU);
            let direction = if direction >= TAU { 0.0 } else { direction };
            let start = SurfacePoint::Face {
                face: f as u32,
                bary: [third; 3],
            };
            let trace =
                self.tri
                    .trace(start, direction, d[0].hypot(d[1]), &TraceOptions::default())?;
            match (trace.stop, trace.end) {
                (TraceStop::Completed, SurfacePoint::Face { face, bary }) => (face as usize, bary),
                (TraceStop::Boundary, SurfacePoint::Edge { edge, .. }) => {
                    // The circumcentre is outside the surface: split the boundary edge the
                    // path leaves through, and look at the face again.
                    if self.split_boundary(edge as usize)? {
                        self.retry(f);
                    } else {
                        self.block(
                            f,
                            BelowReason::NotInsertable {
                                detail: format!("boundary edge {edge} could not be split"),
                            },
                        );
                    }
                    return Ok(());
                }
                (stop, end) => {
                    self.block(
                        f,
                        BelowReason::NotInsertable {
                            detail: format!(
                                "the path to the circumcentre ended at {end} ({stop:?})"
                            ),
                        },
                    );
                    return Ok(());
                }
            }
        };
        match self.examine(g, bary, shortest) {
            Verdict::Clear => {}
            Verdict::Encroaches(edge) => {
                if self.split_boundary(edge as usize)? {
                    self.retry(f);
                } else {
                    self.block(
                        f,
                        BelowReason::NotInsertable {
                            detail: format!("boundary edge {edge} could not be split"),
                        },
                    );
                }
                return Ok(());
            }
            Verdict::Needle(vertex) => {
                let angle_sum = self.tri.angle_sums[vertex as usize];
                self.block(f, BelowReason::NeedleVertex { vertex, angle_sum });
                return Ok(());
            }
        }
        match self.tri.insert_vertex(g, bary) {
            Ok(insertion) => {
                self.n_face_insertions += 1;
                self.after_insertion(&insertion);
            }
            Err(IntrinsicError::InsertionRefused {
                reason:
                    InsertRefusal::TooThin {
                        nearest_edge: Some(edge),
                        ..
                    },
                ..
            }) => {
                // The circumcentre lies on an edge of the face it is in: split that edge at
                // the point's foot (a boundary edge: by the boundary rule).
                let edge = edge as usize;
                if self.tri.is_boundary_edge(edge) {
                    if self.split_boundary(edge)? {
                        self.retry(f);
                        return Ok(());
                    }
                } else {
                    let k = (0..3)
                        .find(|&k| self.tri.halfedge_edge[3 * g + k] as usize == edge)
                        .unwrap_or(0);
                    let from_tail = bary[(k + 1) % 3] / (bary[k] + bary[(k + 1) % 3]);
                    let forward = self.tri.edge_halfedge[edge] as usize == 3 * g + k;
                    let t = if forward { from_tail } else { 1.0 - from_tail };
                    match self.tri.split_edge(edge, t) {
                        Ok(insertion) => {
                            self.n_edge_splits += 1;
                            self.after_insertion(&insertion);
                            return Ok(());
                        }
                        Err(IntrinsicError::InsertionRefused { .. }) => {}
                        Err(other) => return Err(other),
                    }
                }
                self.block(
                    f,
                    BelowReason::NotInsertable {
                        detail: format!(
                            "the circumcentre lies on edge {edge}, which could not be split there"
                        ),
                    },
                );
            }
            Err(IntrinsicError::InsertionRefused { reason, .. }) => {
                self.block(
                    f,
                    BelowReason::NotInsertable {
                        detail: reason.to_string(),
                    },
                );
            }
            Err(other) => return Err(other),
        }
        Ok(())
    }

    /// Puts face `f` back on the heap if it is unchanged and still needs refinement (a face
    /// that changed was put there when it changed). A face that comes back unchanged too
    /// often, each time after a boundary edge was split for it, is given up.
    fn retry(&mut self, f: usize) {
        const MAX_RETRIES: u32 = 64;
        let stamp = self.stamps[f];
        let entry = self.retries.entry(f as u32).or_insert((stamp, 0));
        if entry.0 != stamp {
            *entry = (stamp, 0);
        }
        entry.1 += 1;
        if entry.1 > MAX_RETRIES {
            self.block(
                f,
                BelowReason::NotInsertable {
                    detail: format!(
                        "its circumcentre is still not insertable after {MAX_RETRIES} splits of \
                         boundary edges"
                    ),
                },
            );
            return;
        }
        if self.needs_refinement(f) {
            self.heap
                .push((self.tri.face_area(f).to_bits(), f as u32, self.stamps[f]));
        }
    }

    fn run(&mut self) -> Result<RefineStatus, IntrinsicError> {
        loop {
            if self.n_inserted() >= self.budget {
                return Ok(RefineStatus::InsertionLimit);
            }
            if let Some(e) = self.encroached.pop() {
                if (e as usize) < self.tri.n_edges() && self.is_encroached(e as usize) {
                    self.split_boundary(e as usize)?;
                }
                continue;
            }
            let Some((_, f, stamp)) = self.heap.pop() else {
                return Ok(RefineStatus::Done);
            };
            let f = f as usize;
            if self.stamps[f] != stamp || !self.needs_refinement(f) {
                continue;
            }
            if self
                .blocked
                .get(&(f as u32))
                .is_some_and(|(s, _)| *s == stamp)
            {
                continue;
            }
            self.refine_face(f)?;
        }
    }
}

impl IntrinsicTriangulation {
    /// Delaunay refinement (paper Sec. 4.2; see the [module documentation](self) and SPEC N1,
    /// N4 for the rules): flips to Delaunay, then inserts circumcentres and splits boundary
    /// edges until no triangle has a corner below `options.min_angle` (or a circumradius above
    /// `options.max_circumradius`), except for the corners listed in the report.
    ///
    /// The result is a Delaunay triangulation of the same surface with more vertices. The
    /// new vertices are appended; their places on the input mesh are
    /// [`vertex_location`](Self::vertex_location).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`] for an angle bound that is not in `[0, 35°]`, a
    /// circumradius bound that is not positive, or an invalid Delaunay tolerance. Reaching the
    /// insertion limit is not an error.
    pub fn refine(&mut self, options: &RefineOptions) -> Result<RefineReport, IntrinsicError> {
        let invalid = |reason: String| Err(IntrinsicError::InvalidOptions { reason });
        if !(options.min_angle >= 0.0 && options.min_angle <= 35f64.to_radians()) {
            return invalid(format!(
                "the angle bound of refinement must be between 0 and 35 degrees, got {} degrees",
                options.min_angle.to_degrees()
            ));
        }
        if let Some(r) = options.max_circumradius
            && !(r.is_finite() && r > 0.0)
        {
            return invalid(format!(
                "the circumradius bound of refinement must be positive and finite, got {r}"
            ));
        }
        let first = self.flip_to_delaunay(&crate::DelaunayOptions {
            criterion: options.criterion,
            ..Default::default()
        })?;
        let budget = match options.max_insertions {
            InsertionLimit::PerVertex(k) => k.saturating_mul(self.n_vertices()).max(1000),
            InsertionLimit::Vertices(n) => n,
            InsertionLimit::Unlimited => usize::MAX,
        };
        // Boundary vertices with a small angle between their two boundary edges: the unit of
        // their split distances is the shorter of the two.
        let mut shells = HashMap::new();
        for v in 0..self.n_vertices() {
            if !self.is_boundary_vertex(v) || self.angle_sums[v] >= FRAC_PI_3 {
                continue;
            }
            let first_out = self.vertex_halfedge[v] as usize;
            let last_out = self.outgoing_halfedges(v).last().unwrap_or(first_out);
            let arriving = self.prev(last_out);
            let unit = self.edge_lengths[self.halfedge_edge[first_out] as usize]
                .min(self.edge_lengths[self.halfedge_edge[arriving] as usize]);
            shells.insert(v as u32, unit);
        }
        let n_faces = self.n_faces();
        let n_edges = self.n_edges();
        let mut refiner = Refiner {
            tri: self,
            options: *options,
            queue: EdgeQueue::empty(n_edges),
            heap: BinaryHeap::new(),
            stamps: vec![0; n_faces],
            encroached: Vec::new(),
            blocked: HashMap::new(),
            shells,
            touched: Vec::new(),
            mark: Vec::new(),
            search: 0,
            retries: HashMap::new(),
            budget,
            n_face_insertions: 0,
            n_edge_splits: 0,
            n_boundary_splits: 0,
            n_flips: first.n_flips,
            n_snapped: 0,
        };
        for f in 0..n_faces {
            refiner.touch(f);
        }
        let status = refiner.run()?;

        // What is left.
        let mut below_bound = Vec::new();
        let mut too_large = Vec::new();
        let mut min_angle = f64::INFINITY;
        let mut max_angle = 0.0f64;
        let mut min_angle_reachable = f64::INFINITY;
        for f in 0..refiner.tri.n_faces() {
            if refiner.too_large(f) {
                too_large.push(f as u32);
            }
            let blocked = refiner
                .blocked
                .get(&(f as u32))
                .filter(|(stamp, _)| *stamp == refiner.stamps[f])
                .map(|(_, reason)| reason.clone());
            for k in 0..3 {
                let vertex = refiner.tri.faces[f][k];
                let angle = refiner.tri.corner_angle(3 * f + k);
                min_angle = min_angle.min(angle);
                max_angle = max_angle.max(angle);
                if !refiner.exempt(vertex as usize) {
                    min_angle_reachable = min_angle_reachable.min(angle);
                }
                if angle >= options.min_angle {
                    continue;
                }
                let reason = if refiner.exempt(vertex as usize) {
                    BelowReason::AngleSum {
                        angle_sum: refiner.tri.angle_sums[vertex as usize],
                    }
                } else if let Some(reason) = &blocked {
                    reason.clone()
                } else {
                    BelowReason::InsertionLimit
                };
                below_bound.push(BelowBound {
                    face: f as u32,
                    vertex,
                    angle,
                    reason,
                });
            }
        }
        let report = RefineReport {
            status,
            n_inserted: refiner.n_inserted(),
            n_face_insertions: refiner.n_face_insertions,
            n_edge_splits: refiner.n_edge_splits,
            n_boundary_splits: refiner.n_boundary_splits,
            n_flips: refiner.n_flips,
            n_snapped: refiner.n_snapped,
            min_angle,
            max_angle,
            min_angle_reachable,
            below_bound,
            too_large,
            n_not_delaunay: (0..refiner.tri.n_edges())
                .filter(|&e| refiner.tri.fails_delaunay(e, &options.criterion))
                .count(),
            n_encroached: (0..refiner.tri.n_edges())
                .filter(|&e| refiner.is_encroached(e))
                .count(),
        };
        Ok(report)
    }
}
