//! Moving an inserted vertex inside its flat star without changing the connectivity (SPEC
//! §4.11, I29): the signposts' local update (Sharp, Soliman, Crane 2019, Sec. 3.3.3) on integer
//! coordinates, which the optimal Delaunay iteration uses.
//!
//! The star of an inserted interior vertex `v` is flat (angle sum `2π`), so it unfolds into a
//! polygon `P` around `v` (directions from [`AngleFrames`], SPEC I9). Moving `v` to a point `v'`
//! of the kernel of `P` keeps every triangle of the star; only the spokes change. The input edges
//! that cross `P` are straight chords of `P` between points of its boundary (the link): crossings
//! of link edges, or link vertices that emanate arcs; they are read off the integers by joining
//! the arcs of the star's faces ([`crate::insert`]) across the spokes, plus, when `v` lies on an
//! input edge (two spokes with `n = −1`), that input edge as the chord between the two far ends
//! of those spokes. Then, **exactly from the integers and one side test per chord**:
//!
//! * `n` of the new spoke `v' → u` is the number of chords that separate `v'` from the link
//!   vertex `u` (a chord that ends at `u` does not): the side of `u` is read off the order of the
//!   chord's ends on the link, the side of `v'` is the sign of one orientation test; a `v'`
//!   within `1e-12` mean spoke lengths of a chord is refused ([`IntrinsicError::InsertionRefused`],
//!   the iteration halves the step), so every decision is clear;
//! * the roundabout of `u → v'` for `u ∈ V⁰` follows from that of the link halfedge clockwise of
//!   it by `ic-roundabout-update` with the emanating count of the new face (Eq. 2), as in a flip;
//! * the lengths of the spokes are the signposts' (`it-moved-length`, the law of cosines in the
//!   unfolded star), the refusals those of [`check_star_move`] (turned or too thin triangles);
//! * `q⁰` of `v'` is the point query of the old state at the end of the path from `v`
//!   ([`IntegerTriangulation::point_on_input`], App. A); a place on an input edge or vertex is
//!   refused.
//!
//! Every new face is checked for Eq. 1's parity before anything changes. Unlike
//! [`IntegerTriangulation::move_vertex`] (Sec. 3.6: trace, insert, remove) this never needs a
//! removal, so it is not refused by the straight stars of I27 (every vertex an edge split made,
//! which is what the optimal Delaunay iteration's splits make).

use std::collections::BTreeMap as HashMap;
use std::f64::consts::TAU;

use ddg_intrinsic::formulas::moved_length;
use ddg_intrinsic::reposition::{VertexMove, check_star_move};
use ddg_intrinsic::{
    AngleFrames, InsertRefusal, IntrinsicError, SurfacePoint, TraceOptions, TraceStop,
};

use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{face_emanating, face_valid, n_minus, roundabout_update};
use crate::insert::{ExtractionCache, FaceArcs, RegionPoint};
use crate::recover::recover_in_input_face;

/// A point of the link of the star where a chord ends: link vertex `u_edge` (`sub = 0`) or
/// crossing `sub − 1` of link edge `edge` (from `u_edge` towards `u_{edge+1}`). The order of the
/// keys is the counter-clockwise order on the link.
type LinkKey = (usize, usize);

/// An end of a piece of a chord inside one face of the star.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    Link(LinkKey),
    /// Crossing `index` (counted from `v`) of spoke `spoke`.
    Spoke {
        spoke: usize,
        index: usize,
    },
    /// The vertex `v` itself (an input edge through `v` that emanates from it into a face).
    Center,
}

/// A chord of the star polygon: its two ends on the link, their positions in the layout, and
/// the input edge it is a piece of.
struct StarChord {
    ends: [LinkKey; 2],
    at: [[f64; 2]; 2],
    input_edge: u32,
}

/// Whether `x` lies strictly between `a` and `b` going counter-clockwise from `a`.
fn strictly_between(a: LinkKey, b: LinkKey, x: LinkKey) -> bool {
    if a < b {
        a < x && x < b
    } else {
        x > a || x < b
    }
}

impl IntegerTriangulation {
    /// Moves the inserted interior vertex `v` by `displacement` (in its [`AngleFrames`] frame)
    /// inside its star, keeping the connectivity (see the [module documentation](self)). The
    /// refusals of [`check_star_move`] with `min_shape`, and a new place within rounding of an
    /// input edge, change nothing.
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::MoveRefused`] (input vertex, self-edge, not finite, leaves the star,
    /// too thin), [`IntrinsicError::InsertionRefused`] (the new place within rounding of an input
    /// edge or vertex), [`IntrinsicError::MoveNotLocated`] (the path from `v` stops early);
    /// [`IntegerError::Unsupported`] for a vertex on the boundary (SPEC I7); the errors of
    /// [`point_on_input`](Self::point_on_input) and of the arcs of the star's faces;
    /// [`IntegerError::InvalidCoordinates`] if the integers of the star are not those of a
    /// realization. Nothing is changed in every case.
    pub fn move_vertex_in_star(
        &mut self,
        v: usize,
        displacement: [f64; 2],
        min_shape: f64,
    ) -> Result<VertexMove, IntegerError> {
        let nv = self.complex.n_vertices();
        if v >= nv {
            return Err(IntegerError::IndexOutOfRange {
                what: "vertex",
                index: v,
                bound: nv,
            });
        }
        if !self.is_input_vertex(v) && self.complex.is_boundary_vertex(v) {
            return Err(IntegerError::Unsupported {
                operation: "move_vertex_in_star",
                reason: format!(
                    "vertex {v} is an inserted vertex on the boundary; those are not moved (SPEC \
                     I7)"
                ),
            });
        }
        let c = &self.complex;
        let o: Vec<usize> = c.outgoing_halfedges(v).collect();
        let d = o.len();
        let frames = AngleFrames::of(self);
        let phi = |h: usize| frames.signpost(h);
        let star = check_star_move(self, v, displacement, min_shape, &phi)
            .map_err(IntegerError::Intrinsic)?;
        let from = self.locations[v];
        let r = displacement[0].hypot(displacement[1]);
        let unchanged = VertexMove {
            vertex: v as u32,
            from,
            to: from,
            distance: 0.0,
            back_direction: 0.0,
            snapped: None,
            within_height_bound: true,
            min_shape: star.min_shape,
        };
        if r == 0.0 {
            return Ok(unchanged);
        }
        let phi_v = displacement[1].atan2(displacement[0]);
        // The unfolded star: spoke j to u_j at a[j]; face j = (v, u_j, u_{j+1}) with u_{j+1}
        // at b[j] in its own layout.
        let len = |h: usize| self.lengths[c.edge(h)];
        let polar = |l: f64, angle: f64| [l * angle.cos(), l * angle.sin()];
        let a: Vec<[f64; 2]> = o.iter().map(|&h| polar(len(h), phi(h))).collect();
        let b: Vec<[f64; 2]> = o
            .iter()
            .map(|&h| polar(len(c.prev(h)), phi(h) + self.corner_angle_of(h)))
            .collect();
        let new_lengths: Vec<f64> = o
            .iter()
            .map(|&h| moved_length(len(h), r, phi(h) - phi_v))
            .collect();
        let mean_spoke = o.iter().map(|&h| len(h)).sum::<f64>() / d as f64;

        // The new place: the end of the path from v, located on the input mesh.
        let direction = phi_v.rem_euclid(TAU);
        let direction = if direction >= TAU { 0.0 } else { direction };
        let trace = frames
            .trace(
                SurfacePoint::Vertex(v as u32),
                direction,
                r,
                &TraceOptions::default(),
            )
            .map_err(IntegerError::Intrinsic)?;
        drop(frames);
        if trace.stop != TraceStop::Completed {
            return Err(IntegerError::Intrinsic(IntrinsicError::MoveNotLocated {
                vertex: v as u32,
                detail: format!(
                    "the path from vertex {v} stopped ({:?}) at {} after {:?} of {r:?}",
                    trace.stop, trace.end, trace.length
                ),
            }));
        }
        let refused = |at: SurfacePoint, detail: String| {
            IntegerError::Intrinsic(IntrinsicError::InsertionRefused {
                at,
                reason: InsertRefusal::NotLocated { detail },
            })
        };
        // The chords of the star polygon.
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Vertex(v as u32), detail)
        };
        let mut cache = ExtractionCache::new();
        let mut arcs_of: HashMap<usize, FaceArcs> = HashMap::new();
        let mut pieces: Vec<[End; 2]> = Vec::new();
        let mut piece_edge: Vec<u32> = Vec::new();
        let mut at_link: HashMap<LinkKey, [f64; 2]> = HashMap::new();
        // No input edge crosses a star whose edges all have n = 0 (it lies in one input face):
        // its arcs are not extracted.
        let quiet = o
            .iter()
            .all(|&h| self.normal[c.edge(h)] == 0 && self.normal[c.edge(c.next(h))] == 0);
        for (j, &h) in o.iter().enumerate() {
            if quiet {
                break;
            }
            let f = h / 3;
            let k0 = h % 3;
            let arcs = arcs_of
                .entry(f)
                .or_insert(self.face_arcs_cached(f, &mut cache)?);
            let next = (j + 1) % d;
            let n_next = arcs.n((k0 + 2) % 3);
            at_link.insert((j, 0), a[j]);
            let link = &arcs.sides[(k0 + 1) % 3];
            for (p, &s) in link.params.iter().enumerate() {
                at_link.insert(
                    (j, p + 1),
                    [
                        a[j][0] + s * (b[j][0] - a[j][0]),
                        a[j][1] + s * (b[j][1] - a[j][1]),
                    ],
                );
            }
            let end_of = |rp: RegionPoint| -> Result<End, IntegerError> {
                Ok(match rp {
                    RegionPoint::Corner(k) if k == (k0 + 1) % 3 => End::Link((j, 0)),
                    RegionPoint::Corner(k) if k == (k0 + 2) % 3 => End::Link((next, 0)),
                    RegionPoint::Corner(_) => End::Center,
                    RegionPoint::Crossing(k, p) if k == k0 => End::Spoke { spoke: j, index: p },
                    RegionPoint::Crossing(k, p) if k == (k0 + 1) % 3 => End::Link((j, p + 1)),
                    RegionPoint::Crossing(_, p) => End::Spoke {
                        spoke: next,
                        index: n_next - 1 - p,
                    },
                })
            };
            for chord in arcs.chords() {
                let (x, y) = arcs.chord_ends(chord);
                // Every arc ends at a crossing at least once: its input edge.
                let input_edge = [x, y]
                    .into_iter()
                    .find_map(|rp| match rp {
                        RegionPoint::Crossing(k, p) => match arcs.sides[k].points[p] {
                            SurfacePoint::Edge { edge, .. } => Some(edge),
                            _ => None,
                        },
                        RegionPoint::Corner(_) => None,
                    })
                    .ok_or_else(|| invalid(format!("an arc of face {f} has no crossing")))?;
                pieces.push([end_of(x)?, end_of(y)?]);
                piece_edge.push(input_edge);
            }
        }
        let mut at_spoke: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (i, piece) in pieces.iter().enumerate() {
            for e in piece {
                if let End::Spoke { spoke, index } = *e {
                    at_spoke.entry((spoke, index)).or_default().push(i);
                }
            }
        }
        let mut used = vec![false; pieces.len()];
        let mut chords: Vec<StarChord> = Vec::new();
        // Halves of the input edge through v: from a link end to v.
        let mut through_v: Vec<(LinkKey, u32)> = Vec::new();
        for start in 0..pieces.len() {
            if used[start] {
                continue;
            }
            let Some(side) = pieces[start].iter().position(|e| matches!(e, End::Link(_))) else {
                continue;
            };
            let End::Link(first) = pieces[start][side] else {
                continue;
            };
            let (mut piece, mut out) = (start, pieces[start][1 - side]);
            used[piece] = true;
            let last = loop {
                match out {
                    End::Link(key) => break Some(key),
                    End::Center => break None,
                    End::Spoke { spoke, index } => {
                        let next = at_spoke
                            .get(&(spoke, index))
                            .and_then(|list| list.iter().copied().find(|&q| q != piece))
                            .ok_or_else(|| {
                                invalid(format!(
                                    "crossing {index} of spoke {spoke} ends one arc only"
                                ))
                            })?;
                        if used[next] {
                            return Err(invalid("an input edge closes up inside the star".into()));
                        }
                        used[next] = true;
                        let k = usize::from(pieces[next][0] == out);
                        out = pieces[next][k];
                        piece = next;
                    }
                }
            };
            match last {
                Some(last) => chords.push(StarChord {
                    ends: [first, last],
                    at: [at_link[&first], at_link[&last]],
                    input_edge: piece_edge[start],
                }),
                None => through_v.push((first, piece_edge[start])),
            }
        }
        if let Some(i) = used.iter().position(|u| !u) {
            return Err(invalid(format!(
                "the arc {i} of the star crosses spokes only (a closed input curve around {v})"
            )));
        }
        // An input edge through v: along spokes with n = −1, or emanating from v into faces.
        let mut halves = through_v;
        for j in (0..d).filter(|&j| self.normal[c.edge(o[j])] < 0) {
            let SurfacePoint::Edge { edge, .. } = from else {
                return Err(invalid(format!(
                    "spoke {j} of vertex {v} runs along an input edge, but the vertex is at {from}"
                )));
            };
            halves.push(((j, 0), edge));
        }
        match halves[..] {
            [] => {}
            [(x, g), (y, h)] if g == h => chords.push(StarChord {
                ends: [x, y],
                at: [at_link[&x], at_link[&y]],
                input_edge: g,
            }),
            _ => {
                return Err(invalid(format!(
                    "{} input edge halves end at vertex {v} ({halves:?})",
                    halves.len()
                )));
            }
        }

        let end = trace.end;
        let to = match end {
            SurfacePoint::Vertex(u) => {
                return Err(refused(end, format!("the new place is vertex {u}")));
            }
            // The point query of a point inside a face (SPEC §4.9), with the arcs of the star
            // already extracted.
            SurfacePoint::Face { face, bary }
                if bary.iter().all(|&x| x > 0.0) && o.iter().any(|&h| h / 3 == face as usize) =>
            {
                let f = face as usize;
                if let std::collections::btree_map::Entry::Vacant(slot) = arcs_of.entry(f) {
                    slot.insert(self.face_arcs_cached(f, &mut cache)?);
                }
                let arcs = &arcs_of[&f];
                let region = arcs.choose_region(bary).region;
                recover_in_input_face(self, arcs, region, bary)?.on_input(self, face as usize)?
            }
            _ => self.point_on_input(end)?,
        };
        // On an input edge (the point query snaps within rounding, as the signposts snap a
        // moved vertex, SPEC N24): that edge's chord passes through the new place.
        let on_edge = match to {
            SurfacePoint::Face { .. } => None,
            SurfacePoint::Edge { edge, .. } => Some(edge),
            SurfacePoint::Vertex(u) => {
                return Err(refused(end, format!("the new place is input vertex {u}")));
            }
        };

        // The new normal coordinates of the spokes.
        let p = displacement;
        let mut counts = vec![0i64; d];
        let mut through: Option<&StarChord> = None;
        for chord in &chords {
            if Some(chord.input_edge) == on_edge {
                if through.is_some() {
                    return Err(refused(
                        end,
                        format!("input edge {} crosses the star twice", chord.input_edge),
                    ));
                }
                through = Some(chord);
                continue;
            }
            let [x, y] = chord.at;
            let e = [y[0] - x[0], y[1] - x[1]];
            let side = e[0] * (p[1] - x[1]) - e[1] * (p[0] - x[0]);
            let length = e[0].hypot(e[1]);
            if side.is_nan() || side.abs() <= 1e-12 * mean_spoke * length {
                return Err(refused(
                    end,
                    format!(
                        "the new place is within rounding of the input edge from {:?} to {:?} of \
                         the star",
                        chord.ends[0], chord.ends[1]
                    ),
                ));
            }
            let left = side > 0.0;
            for (i, count) in counts.iter_mut().enumerate() {
                let key = (i, 0);
                if chord.ends.contains(&key) {
                    continue;
                }
                // The link counter-clockwise from the chord's first end to its second lies on
                // the chord's right.
                let u_left = !strictly_between(chord.ends[0], chord.ends[1], key);
                if u_left != left {
                    *count += 1;
                }
            }
        }
        if let Some(g) = on_edge {
            let Some(chord) = through else {
                return Err(refused(
                    end,
                    format!("the new place is on input edge {g}, which does not cross the star"),
                ));
            };
            // The spokes to the ends of that chord at link vertices run along it.
            for (i, count) in counts.iter_mut().enumerate() {
                if chord.ends.contains(&(i, 0)) {
                    *count = -1;
                }
            }
        }
        // The new faces, their parity, and the roundabouts of the halfedges u_j → v.
        let mut roundabouts = vec![NO_ROUNDABOUT; d];
        for (j, &h) in o.iter().enumerate() {
            let k0 = h % 3;
            let mut n = [0i64; 3];
            n[k0] = counts[j];
            n[(k0 + 1) % 3] = self.normal[c.edge(c.next(h))];
            n[(k0 + 2) % 3] = counts[(j + 1) % d];
            if !face_valid(n) {
                return Err(invalid(format!(
                    "the face {} would get the side counts {n:?}, which are not valid",
                    h / 3
                )));
            }
            let u = c.to_vertex(h);
            if self.is_input_vertex(u) {
                let clockwise = c.next(h);
                let deg = self.local.degree(u);
                let rc = self.roundabouts[clockwise];
                let x = (rc < deg)
                    .then(|| {
                        roundabout_update(
                            i64::from(rc),
                            face_emanating(n, (k0 + 1) % 3),
                            n_minus(self.normal[c.edge(clockwise)]),
                            self.local.modulus(u),
                        )
                    })
                    .flatten();
                roundabouts[j] = match x {
                    Some(x) if (0..i64::from(deg)).contains(&x) => x as u32,
                    _ => {
                        return Err(invalid(format!(
                            "the roundabout {rc} of halfedge {clockwise} at input vertex {u} gives \
                             no roundabout for the moved spoke"
                        )));
                    }
                };
            }
        }

        // Change the state.
        for (j, &h) in o.iter().enumerate() {
            let e = self.complex.edge(h);
            self.lengths[e] = new_lengths[j];
            self.normal[e] = counts[j];
            if roundabouts[j] != NO_ROUNDABOUT {
                let t = self
                    .complex
                    .twin(h)
                    .expect("an interior vertex has no boundary spoke");
                self.roundabouts[t] = roundabouts[j];
            }
        }
        self.locations[v] = to;
        Ok(VertexMove {
            to,
            distance: r,
            back_direction: f64::NAN,
            within_height_bound: star.within_height_bound,
            ..unchanged
        })
    }

    /// The corner angle at the tail of halfedge `h`, from the lengths.
    fn corner_angle_of(&self, h: usize) -> f64 {
        ddg_intrinsic::Triangulation::corner_angle(self, h)
    }
}
