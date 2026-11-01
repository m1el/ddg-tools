//! Straight segments of the triangulation as paths over the input mesh (SPEC §4.10, I29): the
//! way back of `ddg_intrinsic`'s simplicial export, with the signature, output and errors of
//! [`ddg_intrinsic::IntrinsicTriangulation::segments_on_input`].
//!
//! A straight segment inside a face `f` of `T¹` crosses the input edges that cross `f`. Those
//! are the arcs of the face ([`crate::insert`]): chords of the face laid out in the plane,
//! between the crossings of its sides (positions from `ExtractIntrinsicEdge`) or from an
//! emanating corner, given by the integers alone. The arcs cut the face into regions, each
//! inside one input face (the region's input face, as for App. A, [`crate::recover`]); the
//! regions and arcs form a tree (every arc is a chord of a convex triangle, so it separates the
//! face in two). **Which input edges the segment crosses, and in which order, is decided
//! combinatorially:** the region of each end (an end inside the face: the region the face split
//! would pick, I15; on a side: the gap between the side's crossings it lies in, counted as the
//! point query counts it; at a corner: the region at the corner towards the other end), then the
//! arcs on the tree path between the two regions, in order. Floating point decides only where
//! the segment meets each arc (the intersection of two lines in the face's layout, clamped into
//! the arc, I5) and so the parameter on the input edge (interpolated between the arc's ends).
//!
//! A segment along a side of the face is the piece of that side's path over the input mesh
//! between its ends (`ExtractIntrinsicEdge`), or of the input edge the side runs along.
//!
//! **Checks** (the certificate's spirit, SPEC §8): every piece's two end points must lie in the
//! closure of the piece's input face, and the length of the path measured piece by piece in the
//! input faces must equal the length of the segment within `options.tolerance` times that length
//! ([`PathOnInput::mismatch`] is the difference, an absolute length). Either failure is
//! [`IntrinsicError::SegmentTrace`], as a trace that does not arrive is for the signposts.

use std::collections::{HashMap, VecDeque};

use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic::{
    EdgeTraceOptions, IntrinsicError, PathOnInput, PathsOnInput, Segment, SurfacePoint,
};

use crate::coords::IntegerTriangulation;
use crate::error::IntegerError;
use crate::insert::{Chord, FaceArcs, Region, RegionPoint};
use crate::recover::bary_in_input_face;

/// How far the sum of barycentric coordinates may be from one (as `ddg_intrinsic`).
const BARY_SUM_TOLERANCE: f64 = 1e-9;

/// What is wrong with `point` on a surface with the given element counts (the rules of
/// `ddg_intrinsic`'s point checks, so that both representations refuse the same nodes).
fn point_defect(
    point: &SurfacePoint,
    n_vertices: usize,
    n_edges: usize,
    n_faces: usize,
) -> Option<&'static str> {
    match *point {
        SurfacePoint::Vertex(v) => {
            (v as usize >= n_vertices).then_some("the vertex does not exist")
        }
        SurfacePoint::Edge { edge, t } => {
            if edge as usize >= n_edges {
                Some("the edge does not exist")
            } else if !(t.is_finite() && (0.0..=1.0).contains(&t)) {
                Some("the edge parameter is not in [0, 1]")
            } else {
                None
            }
        }
        SurfacePoint::Face { face, bary } => {
            if face as usize >= n_faces {
                Some("the face does not exist")
            } else if !bary.iter().all(|b| b.is_finite() && *b >= 0.0) {
                Some("a barycentric coordinate is negative or not finite")
            } else if ((bary[0] + bary[1] + bary[2]) - 1.0).abs() > BARY_SUM_TOLERANCE {
                Some("the barycentric coordinates do not sum to one")
            } else {
                None
            }
        }
    }
}

/// Where a node of a segment lies in the face of the segment.
#[derive(Debug, Clone, Copy)]
enum NodeKind {
    /// At the corner in this slot.
    Corner(usize),
    /// On the side in this slot, at parameter `t` of its edge (in the edge's direction).
    Side { slot: usize, t: f64 },
    /// Inside the face, with these barycentric coordinates.
    Inside([f64; 3]),
}

impl NodeKind {
    /// The sides of the face the node lies on (a corner: the side leaving it and the side
    /// arriving at it).
    fn sides(self) -> [Option<usize>; 2] {
        match self {
            NodeKind::Corner(k) => [Some(k), Some((k + 2) % 3)],
            NodeKind::Side { slot, .. } => [Some(slot), None],
            NodeKind::Inside(_) => [None, None],
        }
    }
}

/// The arcs of a face with the tree of its regions.
struct FaceTree {
    arcs: FaceArcs,
    /// Per region: the arcs that bound it, each with the region on its other side.
    neighbours: HashMap<Region, Vec<(Chord, Region)>>,
}

impl FaceTree {
    fn new(arcs: FaceArcs) -> FaceTree {
        let mut sides: HashMap<Chord, [Option<Region>; 2]> = HashMap::new();
        for region in arcs.arrangement.regions() {
            for (chord, sign) in arcs.bounds(region) {
                sides.entry(chord).or_default()[usize::from(sign < 0.0)] = Some(region);
            }
        }
        let mut neighbours: HashMap<Region, Vec<(Chord, Region)>> = HashMap::new();
        for (chord, [a, b]) in sides {
            if let (Some(a), Some(b)) = (a, b) {
                neighbours.entry(a).or_default().push((chord, b));
                neighbours.entry(b).or_default().push((chord, a));
            }
        }
        for list in neighbours.values_mut() {
            list.sort_unstable_by_key(|&(chord, _)| match chord {
                Chord::Around(k, t) => (0, k, t),
                Chord::Fan(s) => (1, 0, s),
            });
        }
        FaceTree { arcs, neighbours }
    }

    /// The region that touches the gap `slot` of side `k` (between its crossings `slot − 1`
    /// and `slot` from the side's tail): crossings `[0, c_k)` of side `k` end arcs around corner
    /// `k`, then come the fan arcs if `k` is opposite the emanating corner, then the arcs around
    /// corner `k + 1`.
    fn gap_region(&self, k: usize, slot: usize) -> Region {
        let (ck, ck1) = (self.arcs.c(k), self.arcs.c((k + 1) % 3));
        let n = self.arcs.n(k);
        if slot < ck {
            return Region::Corner {
                corner: k,
                index: slot as i64,
            };
        }
        if slot + ck1 > n {
            return Region::Corner {
                corner: (k + 1) % 3,
                index: (n - slot) as i64,
            };
        }
        self.middle_at(k, slot - ck)
    }

    /// The middle region at offset `s` from the arcs around corner `k` on side `k`.
    fn middle_at(&self, k: usize, s: usize) -> Region {
        let a = &self.arcs.arrangement;
        if a.fan == 0 {
            return Region::Middle { index: 0 };
        }
        let kappa = a.kappa;
        Region::Middle {
            index: if k == (kappa + 1) % 3 {
                s as i64
            } else if k == (kappa + 2) % 3 {
                a.fan
            } else {
                0
            },
        }
    }

    /// The regions that touch corner `k`.
    fn corner_regions(&self, k: usize) -> Vec<Region> {
        let a = &self.arcs.arrangement;
        if self.arcs.c(k) > 0 {
            return vec![Region::Corner {
                corner: k,
                index: 0,
            }];
        }
        if a.fan > 0 && k == a.kappa {
            return (0..=a.fan).map(|s| Region::Middle { index: s }).collect();
        }
        // The middle region at the corner: offset 0 on the side leaving it.
        vec![self.middle_at(k, 0)]
    }

    /// The tree path from `from` to `to`: the arcs crossed and the regions entered, in order.
    fn path(&self, from: Region, to: Region) -> Option<Vec<(Chord, Region)>> {
        let mut previous: HashMap<Region, (Chord, Region)> = HashMap::new();
        let mut queue = VecDeque::from([from]);
        let mut seen = std::collections::HashSet::from([from]);
        while let Some(r) = queue.pop_front() {
            if r == to {
                break;
            }
            for &(chord, next) in self.neighbours.get(&r).map_or(&[][..], |v| &v[..]) {
                if seen.insert(next) {
                    previous.insert(next, (chord, r));
                    queue.push_back(next);
                }
            }
        }
        if !seen.contains(&to) {
            return None;
        }
        let mut out = Vec::new();
        let mut at = to;
        while at != from {
            let (chord, back) = previous[&at];
            out.push((chord, at));
            at = back;
        }
        out.reverse();
        Some(out)
    }

    /// The tree distance from `from` to every region.
    fn distances(&self, from: Region) -> HashMap<Region, usize> {
        let mut out = HashMap::from([(from, 0usize)]);
        let mut queue = VecDeque::from([from]);
        while let Some(r) = queue.pop_front() {
            let d = out[&r];
            for &(_, next) in self.neighbours.get(&r).map_or(&[][..], |v| &v[..]) {
                if let std::collections::hash_map::Entry::Vacant(e) = out.entry(next) {
                    e.insert(d + 1);
                    queue.push_back(next);
                }
            }
        }
        out
    }

    /// The position in the face's layout of a region point.
    fn position(&self, p: RegionPoint) -> [f64; 2] {
        let b = self.arcs.bary_of(p);
        let v = &self.arcs.layout;
        [
            b[0] * v[0][0] + b[1] * v[1][0] + b[2] * v[2][0],
            b[0] * v[0][1] + b[1] * v[1][1] + b[2] * v[2][1],
        ]
    }
}

/// `a × b`.
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

impl IntegerTriangulation {
    /// Carries straight segments of the triangulation to the input mesh (SPEC §4.10, I29): the
    /// contract of [`ddg_intrinsic::IntrinsicTriangulation::segments_on_input`] (the same
    /// nodes, segments, output and errors) on integer coordinates. Each node is placed once
    /// ([`point_on_input`](Self::point_on_input)); the input edges a segment crosses come from
    /// the integers, their positions from the face's layout (see the
    /// [module documentation](crate::segments)). [`PathOnInput::mismatch`] is the difference
    /// between the length of the path measured in the input faces and the segment's length;
    /// above `options.tolerance` times that length it is an error. `options.max_crossings` is
    /// not used (the work is bounded by the integers).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`] for a node that does not exist or does not lie in the
    /// face of a segment that uses it (and for the end of an edge or a corner not given as the
    /// vertex), [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::InvalidOptions`],
    /// [`IntrinsicError::InvalidLength`] for a segment of length zero,
    /// [`IntrinsicError::SegmentTrace`] for a path whose pieces do not lie in their input
    /// faces or whose length differs from the segment's by more than the tolerance; this
    /// crate's errors of the point queries and of `ExtractIntrinsicEdge` (as
    /// [`IntegerError`]).
    pub fn segments_on_input(
        &self,
        nodes: &[SurfacePoint],
        segments: &[Segment],
        options: &EdgeTraceOptions,
    ) -> Result<PathsOnInput, IntegerError> {
        let shared = |e: IntrinsicError| IntegerError::Intrinsic(e);
        if !(options.tolerance.is_finite() && options.tolerance >= 0.0) {
            return Err(shared(IntrinsicError::InvalidOptions {
                reason: format!(
                    "the tolerance of a segment trace must be finite and not negative, got {}",
                    options.tolerance
                ),
            }));
        }
        let c = &self.complex;
        let mut located = Vec::with_capacity(nodes.len());
        for &node in nodes {
            if let Some(reason) = point_defect(&node, c.n_vertices(), c.n_edges(), c.n_faces()) {
                return Err(shared(IntrinsicError::InvalidPoint {
                    point: node,
                    reason,
                }));
            }
            match node {
                SurfacePoint::Edge { t, .. } if t == 0.0 || t == 1.0 => {
                    return Err(shared(IntrinsicError::InvalidPoint {
                        point: node,
                        reason: "the end of an edge must be given as the vertex",
                    }));
                }
                SurfacePoint::Face { bary, .. } if bary.contains(&1.0) => {
                    return Err(shared(IntrinsicError::InvalidPoint {
                        point: node,
                        reason: "the corner of a face must be given as the vertex",
                    }));
                }
                _ => {}
            }
            located.push(self.point_on_input(node)?);
        }
        let mut trees: HashMap<usize, FaceTree> = HashMap::new();
        let mut paths = Vec::with_capacity(segments.len());
        let mut max_relative_mismatch = 0.0f64;
        for segment in segments {
            let f = segment.face as usize;
            if f >= c.n_faces() {
                return Err(shared(IntrinsicError::IndexOutOfRange {
                    what: "face",
                    index: f,
                    bound: c.n_faces(),
                }));
            }
            for &n in &segment.nodes {
                if n as usize >= nodes.len() {
                    return Err(shared(IntrinsicError::IndexOutOfRange {
                        what: "node",
                        index: n as usize,
                        bound: nodes.len(),
                    }));
                }
            }
            let [from, to] = segment.nodes.map(|n| n as usize);
            let l = self.face_lengths(f);
            let v = triangle_layout(l).ok_or(shared(IntrinsicError::DegenerateFace {
                surface: "intrinsic triangulation",
                face: f as u32,
                lengths: l,
            }))?;
            let (p, kind) = self.node_in_face(f, &v, nodes[from])?;
            let (q, kind_to) = self.node_in_face(f, &v, nodes[to])?;
            let d = [q[0] - p[0], q[1] - p[1]];
            let length = d[0].hypot(d[1]);
            if length <= 0.0 || length.is_nan() {
                return Err(shared(IntrinsicError::InvalidLength { length }));
            }
            let tolerance = options.tolerance * length;
            let ends = ((kind, p, located[from]), (kind_to, q, located[to]));
            let shared_side = kind
                .sides()
                .into_iter()
                .flatten()
                .find(|s| kind_to.sides().contains(&Some(*s)));
            let (points, faces) = match shared_side {
                Some(side) => self.path_along_side(f, side, ends)?,
                None => {
                    if let std::collections::hash_map::Entry::Vacant(e) = trees.entry(f) {
                        e.insert(FaceTree::new(self.face_arcs(f)?));
                    }
                    self.path_across(f, &trees[&f], ends)?
                }
            };
            let mismatch = match self.path_length(&points, &faces) {
                Some(total) => (total - length).abs(),
                None => {
                    return Err(shared(IntrinsicError::SegmentTrace {
                        face: f as u32,
                        length,
                        tolerance,
                        mismatch: None,
                    }));
                }
            };
            if !(mismatch <= tolerance) {
                return Err(shared(IntrinsicError::SegmentTrace {
                    face: f as u32,
                    length,
                    tolerance,
                    mismatch: Some(mismatch),
                }));
            }
            max_relative_mismatch = max_relative_mismatch.max(mismatch / length);
            paths.push(PathOnInput {
                points,
                faces,
                mismatch,
            });
        }
        Ok(PathsOnInput {
            nodes: located,
            paths,
            max_relative_mismatch,
        })
    }

    /// The position of a node in the layout `v` of face `f` and where it lies in the face (as
    /// the signposts' `node_in_face`: a vertex at the first corner that is it, a point of an
    /// edge on the first side that is it).
    fn node_in_face(
        &self,
        f: usize,
        v: &[[f64; 2]; 3],
        node: SurfacePoint,
    ) -> Result<([f64; 2], NodeKind), IntegerError> {
        let c = &self.complex;
        let outside = IntegerError::Intrinsic(IntrinsicError::InvalidPoint {
            point: node,
            reason: "the point is not a point of the face of its segment",
        });
        match node {
            SurfacePoint::Vertex(u) => {
                let k = (0..3).find(|&k| c.faces()[f][k] == u).ok_or(outside)?;
                Ok((v[k], NodeKind::Corner(k)))
            }
            SurfacePoint::Edge { edge, t } => {
                let k = (0..3)
                    .find(|&k| c.edge(3 * f + k) == edge as usize)
                    .ok_or(outside)?;
                let forward = c.edge_halfedge(edge as usize) == 3 * f + k;
                let th = if forward { t } else { 1.0 - t };
                let (a, b) = (v[k], v[(k + 1) % 3]);
                Ok((
                    [a[0] + th * (b[0] - a[0]), a[1] + th * (b[1] - a[1])],
                    NodeKind::Side { slot: k, t },
                ))
            }
            SurfacePoint::Face { face, bary } => {
                if face as usize != f {
                    return Err(outside);
                }
                Ok((
                    [
                        bary[1] * v[1][0] + bary[2] * v[2][0],
                        bary[1] * v[1][1] + bary[2] * v[2][1],
                    ],
                    NodeKind::Inside(bary),
                ))
            }
        }
    }

    /// The path of a segment along side `side` of face `f`: the piece of the input edge the
    /// side runs along, or of the side's own path over the input mesh between the two ends.
    #[allow(clippy::type_complexity)]
    fn path_along_side(
        &self,
        f: usize,
        side: usize,
        ((kind, _, from), (kind_to, _, to)): (
            (NodeKind, [f64; 2], SurfacePoint),
            (NodeKind, [f64; 2], SurfacePoint),
        ),
    ) -> Result<(Vec<SurfacePoint>, Vec<u32>), IntegerError> {
        let c = &self.complex;
        let mesh = self.input.mesh();
        let h = 3 * f + side;
        let e = c.edge(h);
        let forward = c.edge_halfedge(e) == h;
        // The parameter of a node along the edge (in the edge's direction).
        let along = |kind: NodeKind| match kind {
            NodeKind::Corner(k) => {
                let th = if k == side { 0.0 } else { 1.0 };
                if forward { th } else { 1.0 - th }
            }
            NodeKind::Side { t, .. } => t,
            // Not called: a node inside the face lies on no side.
            NodeKind::Inside(_) => 0.0,
        };
        let (tp, tq) = (along(kind), along(kind_to));
        if self.normal[e] < 0 {
            let g = self.input_edge_of_along(e)?;
            return Ok((vec![from, to], vec![(mesh.edge_halfedge(g) / 3) as u32]));
        }
        let trace = self.extract_intrinsic_edge(e, &EdgeTraceOptions::default())?;
        let x = &trace.crossings;
        let len = self.lengths[e];
        let before = |t: f64| x.iter().filter(|c| c.distance < t * len).count();
        let (bp, bq) = (before(tp), before(tq));
        let invalid = |detail: String| {
            IntegerError::invalid(
                crate::error::ValidityItem::Extraction,
                crate::error::Place::Edge(e as u32),
                detail,
            )
        };
        let twin_face = |g: usize| -> Result<u32, IntegerError> {
            mesh.twin(g)
                .map(|t| (t / 3) as u32)
                .ok_or_else(|| invalid(format!("it crosses the boundary input halfedge {g}")))
        };
        // The input face at a point of the edge with `b` crossings before it (as the point
        // query of a point of the edge decides it).
        let face_at = |b: usize| -> Result<u32, IntegerError> {
            if let Some(next) = x.get(b) {
                Ok(next.halfedge / 3)
            } else if let Some(last) = x.last() {
                twin_face(last.halfedge as usize)
            } else {
                trace
                    .first_face
                    .ok_or_else(|| invalid("it has no crossing and no first input face".into()))
            }
        };
        let mut points = vec![from];
        let mut faces = vec![face_at(bp)?];
        if bp <= bq {
            for crossing in &x[bp..bq] {
                points.push(SurfacePoint::Edge {
                    edge: crossing.edge,
                    t: crossing.t,
                });
                faces.push(twin_face(crossing.halfedge as usize)?);
            }
        } else {
            // Backwards along the edge: crossing `i` leads from the twin's face into the face of
            // its halfedge.
            for crossing in x[bq..bp].iter().rev() {
                points.push(SurfacePoint::Edge {
                    edge: crossing.edge,
                    t: crossing.t,
                });
                faces.push(crossing.halfedge / 3);
            }
        }
        points.push(to);
        Ok((points, faces))
    }

    /// The path of a segment of face `f` whose ends do not share a side: the arcs on the tree
    /// path between the regions of its ends.
    #[allow(clippy::type_complexity)]
    fn path_across(
        &self,
        f: usize,
        tree: &FaceTree,
        ((kind, p, from), (kind_to, q, to)): (
            (NodeKind, [f64; 2], SurfacePoint),
            (NodeKind, [f64; 2], SurfacePoint),
        ),
    ) -> Result<(Vec<SurfacePoint>, Vec<u32>), IntegerError> {
        let arcs = &tree.arcs;
        // The regions an end may start in: one, except at a corner that emanates arcs.
        let candidates = |kind: NodeKind| -> Vec<Region> {
            match kind {
                NodeKind::Corner(k) => tree.corner_regions(k),
                NodeKind::Side { slot, t } => {
                    let side = &arcs.sides[slot];
                    let n = arcs.n(slot);
                    let before = match &side.trace {
                        Some(trace) => {
                            let at = t * self.lengths[self.complex.edge(3 * f + slot)];
                            trace.crossings.iter().filter(|c| c.distance < at).count()
                        }
                        None => 0,
                    };
                    let gap = if side.reversed { n - before } else { before };
                    vec![tree.gap_region(slot, gap)]
                }
                NodeKind::Inside(bary) => vec![arcs.choose_region(bary).region],
            }
        };
        let (mut a, mut b) = (candidates(kind), candidates(kind_to));
        // At an emanating corner: the region nearest (in the tree) to the other end's. Two
        // ends at corners share a side, so at most one end has several candidates.
        if a.len() > 1 {
            let d = tree.distances(b[0]);
            a.sort_by_key(|r| d.get(r).copied().unwrap_or(usize::MAX));
        }
        if b.len() > 1 {
            let d = tree.distances(a[0]);
            b.sort_by_key(|r| d.get(r).copied().unwrap_or(usize::MAX));
        }
        let (start, end) = (a[0], b[0]);
        let route = tree.path(start, end).ok_or_else(|| {
            IntegerError::invalid(
                crate::error::ValidityItem::Extraction,
                crate::error::Place::Face(f as u32),
                format!("regions {start} and {end} of face {f} are not connected by its arcs"),
            )
        })?;
        let mut points = vec![from];
        let mut faces = vec![self.region_input_face(arcs, start)? as u32];
        let d = [q[0] - p[0], q[1] - p[1]];
        let mesh = self.input.mesh();
        for &(chord, region) in &route {
            let (ra, rb) = arcs.chord_ends(chord);
            let (pa, pb) = (tree.position(ra), tree.position(rb));
            let e = [pb[0] - pa[0], pb[1] - pa[1]];
            let denom = cross(e, d);
            let u = cross([p[0] - pa[0], p[1] - pa[1]], d) / denom;
            // SPEC I5: clamped into the arc; a segment parallel to it meets it at its middle.
            let u = if u.is_finite() {
                u.clamp(0.0, 1.0)
            } else {
                0.5
            };
            let input_point = |r: RegionPoint| -> Result<SurfacePoint, IntegerError> {
                Ok(match r {
                    RegionPoint::Corner(k) => self.locations[self.complex.faces()[f][k] as usize],
                    RegionPoint::Crossing(k, i) => arcs.sides[k].points[i],
                })
            };
            let (xa, xb) = (input_point(ra)?, input_point(rb)?);
            let SurfacePoint::Edge { edge: g, t: tb } = xb else {
                return Err(self.bad_arc(f, chord, "its end on a side is not a crossing"));
            };
            let ta = match xa {
                SurfacePoint::Edge { edge, t } if edge == g => t,
                SurfacePoint::Vertex(w) => {
                    let [lo, hi] = mesh.edge_vertices(g as usize);
                    if w as usize == lo {
                        0.0
                    } else if w as usize == hi {
                        1.0
                    } else {
                        return Err(self.bad_arc(f, chord, "its corner is not an end of its edge"));
                    }
                }
                _ => return Err(self.bad_arc(f, chord, "its two ends are on different edges")),
            };
            points.push(SurfacePoint::Edge {
                edge: g,
                t: ta + u * (tb - ta),
            });
            faces.push(self.region_input_face(arcs, region)? as u32);
        }
        points.push(to);
        Ok((points, faces))
    }

    fn bad_arc(&self, f: usize, chord: Chord, what: &str) -> IntegerError {
        IntegerError::invalid(
            crate::error::ValidityItem::Extraction,
            crate::error::Place::Face(f as u32),
            format!("the arc {chord:?} of face {f}: {what}"),
        )
    }

    /// The length of a path over the input mesh, measured piece by piece in the layout of each
    /// piece's input face; `None` if an end of a piece is not a point of its face's closure.
    fn path_length(&self, points: &[SurfacePoint], faces: &[u32]) -> Option<f64> {
        let mesh = self.input.mesh();
        let lengths = self.input.geometry().edge_lengths();
        let mut total = 0.0;
        for (i, &g) in faces.iter().enumerate() {
            let g = g as usize;
            let (b0, b1) = (
                bary_in_input_face(self, g, points[i])?,
                bary_in_input_face(self, g, points[i + 1])?,
            );
            let v = triangle_layout([0, 1, 2].map(|s| lengths[mesh.edge(3 * g + s)]))?;
            let at = |b: [f64; 3]| {
                [
                    b[1] * v[1][0] + b[2] * v[2][0],
                    b[1] * v[1][1] + b[2] * v[2][1],
                ]
            };
            let (x0, x1) = (at(b0), at(b1));
            total += (x1[0] - x0[0]).hypot(x1[1] - x0[1]);
        }
        Some(total)
    }
}
