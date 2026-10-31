//! Vertex insertion and edge split (paper Sec. 3.3.2, Alg. 7 and 8).
//!
//! A new vertex is a point of the surface, so it has a place on the input mesh. That place is
//! found before anything is changed: by a tracing query over the input mesh from the nearest
//! corner (Sec. 3.2.3), or, for a split of an edge that lies along an input edge, by
//! interpolating along that input edge (exact; this is always the case on the boundary). The
//! direction back along that trace is the first signpost of the new vertex; the others follow
//! from the corner angles around it.
//!
//! # Frames of inserted vertices (SPEC N10)
//!
//! The signposts of an inserted vertex are angles in the frame of the input element it lies
//! in: the face (angle zero along its first halfedge) or the edge (angle zero along it, from
//! its lower-numbered vertex). An inserted vertex is flat, so no scaling applies. The exception
//! is a vertex on the boundary: its signposts run, as at every boundary vertex, from `0` at
//! the boundary halfedge that leaves it to `π`; where the input boundary edge is numbered
//! against that halfedge the two frames differ by `π`
//! ([`IntrinsicTriangulation::direction_on_input`]).

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::corner_angle;
use ddg_core::topology::NONE;

use crate::error::{InsertRefusal, IntrinsicError};
use crate::flip::triangle_is_valid;
use crate::formulas::{FLIP_SHAPE_LIMIT, corner_distances, signpost_update, triangle_layout};
use crate::point::{SurfacePoint, point_defect};
use crate::surface::Surface;
use crate::trace::{Trace, TraceOptions, TraceStop, run};
use crate::triangulation::IntrinsicTriangulation;

/// A point of an input face counts as lying on an edge of that face when its distance from
/// the edge is at most this fraction of the longest edge of the face (SPEC N10). The place of
/// an inserted vertex is then stored as a point of the input edge. This is the snapping of the
/// paper's Appendix A; it is reported in [`Insertion::snapped`].
///
/// The measure is a distance, not a barycentric coordinate: in a needle of the input mesh
/// (length 0.5, width 1e-4 on one of the test meshes) a point one rounding error of its
/// position away from a long edge has the barycentric coordinate `2e-12`, and a vertex left
/// there makes overlay pieces of the size of a rounding error.
pub const LOCATION_SNAP: f64 = 1e-12;

/// What an insertion or an edge split did.
#[derive(Debug, Clone, PartialEq)]
pub struct Insertion {
    /// The new vertex.
    pub vertex: u32,
    /// Its place on the input mesh.
    pub location: SurfacePoint,
    /// The new edges: from each corner of the face to the vertex (three), or, for a split of
    /// edge `e` from `a` to `b`: the half from the vertex to `b` (the half from `a` keeps the
    /// id `e`), then the edges to the one or two opposite vertices.
    pub new_edges: Vec<u32>,
    /// The new faces (two, or one for a split of a boundary edge). The face the point was in,
    /// or the faces of the split edge, keep their ids.
    pub new_faces: Vec<u32>,
    /// The barycentric coordinate that the traced place had opposite the input edge it was
    /// moved onto ([`LOCATION_SNAP`]); `None` if it was not moved.
    pub snapped: Option<f64>,
}

/// Where a new vertex is on the input mesh and the direction, in the frame of that place,
/// that points from it to the vertex it was traced from.
pub(crate) struct Placed {
    pub(crate) location: SurfacePoint,
    pub(crate) back: f64,
    pub(crate) snapped: Option<f64>,
}

/// What [`IntrinsicTriangulation::snap_onto_input_edge`] made of a point of an input face.
pub(crate) enum EdgeSnap {
    /// Not within the snapping distance of a side (or not a point of a face): unchanged.
    Far,
    /// Moved onto the input edge: its place there, the direction in the frame of the edge,
    /// and the barycentric coordinate it had opposite that edge.
    Onto {
        location: SurfacePoint,
        back: f64,
        coordinate: f64,
    },
    /// Within the snapping distance of a side and exactly at a corner of the face.
    Corner,
    /// The input face has no positive area.
    Degenerate(u32),
}

pub(crate) fn wrap(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    if a >= TAU { 0.0 } else { a }
}

impl IntrinsicTriangulation {
    /// The direction `phi`, given in the frame of vertex `v` (a signpost direction), in the
    /// frame of the place of `v` on the input mesh, as [`InputSurface::trace`] wants it.
    /// They are the same except at an inserted vertex on a boundary edge of the input mesh
    /// whose only halfedge runs from its higher-numbered to its lower-numbered vertex (see the
    /// [module documentation](self)).
    ///
    /// [`InputSurface::trace`]: crate::InputSurface::trace
    pub fn direction_on_input(&self, v: usize, phi: f64) -> f64 {
        match self.locations[v] {
            SurfacePoint::Edge { edge, .. } => {
                let mesh = self.input.mesh();
                let h = mesh.edge_halfedge(edge as usize);
                let backward = mesh.from_vertex(h) != mesh.edge_vertices(edge as usize)[0];
                if mesh.twin(h).is_none() && backward {
                    wrap(phi + PI)
                } else {
                    phi
                }
            }
            _ => phi,
        }
    }

    /// The parameter of vertex `v` along input edge `input_edge`, if it lies on it.
    pub(crate) fn parameter_on_input_edge(&self, v: usize, input_edge: usize) -> Option<f64> {
        match self.locations[v] {
            SurfacePoint::Vertex(u) => {
                let [a, b] = self.input.mesh().edge_vertices(input_edge);
                if u as usize == a {
                    Some(0.0)
                } else if u as usize == b {
                    Some(1.0)
                } else {
                    None
                }
            }
            SurfacePoint::Edge { edge, t } if edge as usize == input_edge => Some(t),
            _ => None,
        }
    }

    /// Traces from vertex `from` in the signpost direction `phi` for `length` over the input
    /// mesh and returns where that is, with the direction back.
    fn place_by_trace(&self, from: usize, phi: f64, length: f64) -> Result<Placed, InsertRefusal> {
        let not_located = |detail: String| InsertRefusal::NotLocated { detail };
        let trace: Trace = run(
            &*self.input,
            self.locations[from],
            self.direction_on_input(from, phi),
            length,
            &TraceOptions::default(),
            None,
        )
        .map_err(|e| not_located(e.to_string()))?
        .trace;
        if trace.stop != TraceStop::Completed {
            return Err(not_located(format!(
                "the path from vertex {from} stopped at {} after a length of {:?} of {length:?}",
                trace.end, trace.length
            )));
        }
        self.place_at_end(&trace, LOCATION_SNAP, f64::INFINITY)
    }

    /// The place on the input mesh at the end of a completed trace over the input mesh, with
    /// the direction back along it; a point within `snap` times the longest edge of its input
    /// face (and within `max_distance`) from an edge of that face ([`LOCATION_SNAP`] for an
    /// insertion, [`MOVE_SNAP`](crate::reposition::MOVE_SNAP) for a move)
    /// face is moved onto that edge (SPEC N10). Refused if the end is an input vertex.
    pub(crate) fn place_at_end(
        &self,
        trace: &Trace,
        snap: f64,
        max_distance: f64,
    ) -> Result<Placed, InsertRefusal> {
        let not_located = |detail: String| InsertRefusal::NotLocated { detail };
        match trace.end {
            SurfacePoint::Vertex(u) => Err(not_located(format!(
                "the point is input vertex {u}, which is a vertex already"
            ))),
            SurfacePoint::Edge { .. } => Ok(Placed {
                location: trace.end,
                back: trace.back_direction,
                snapped: None,
            }),
            SurfacePoint::Face { .. } => {
                match self.snap_onto_input_edge(trace.end, trace.back_direction, snap, max_distance)
                {
                    EdgeSnap::Far => Ok(Placed {
                        location: trace.end,
                        back: trace.back_direction,
                        snapped: None,
                    }),
                    EdgeSnap::Onto {
                        location,
                        back,
                        coordinate,
                    } => Ok(Placed {
                        location,
                        back,
                        snapped: Some(coordinate),
                    }),
                    EdgeSnap::Corner => Err(not_located(format!(
                        "the point {} is a corner of its input face",
                        trace.end
                    ))),
                    EdgeSnap::Degenerate(face) => Err(not_located(format!(
                        "input face {face} has no positive area"
                    ))),
                }
            }
        }
    }

    /// The point `end` of an input face, with the direction `back` in the frame of that face,
    /// as a point of the input edge of the face it is nearest to, if its distance from that
    /// edge is at most `snap` times the longest edge of the face and at most `max_distance`
    /// (SPEC N10, N24, N28). The direction is then given in the frame of the edge. Points
    /// that are not points of a face are [`EdgeSnap::Far`].
    pub(crate) fn snap_onto_input_edge(
        &self,
        end: SurfacePoint,
        back: f64,
        snap: f64,
        max_distance: f64,
    ) -> EdgeSnap {
        let SurfacePoint::Face { face, bary } = end else {
            return EdgeSnap::Far;
        };
        let mesh = self.input.mesh();
        let l = [0, 1, 2]
            .map(|s| self.input.geometry().edge_lengths()[mesh.edge(3 * face as usize + s)]);
        let Some(v) = triangle_layout(l) else {
            return EdgeSnap::Degenerate(face);
        };
        // The distance from the side opposite corner k (halfedge k + 1): the barycentric
        // coordinate times the height over that side.
        let double_area = v[1][0] * v[2][1];
        let distance = [0, 1, 2].map(|k| bary[k] * double_area / l[(k + 1) % 3]);
        let k = (0..3)
            .min_by(|&a, &b| distance[a].total_cmp(&distance[b]))
            .unwrap_or(0);
        let longest = l[0].max(l[1]).max(l[2]);
        if distance[k] > (snap * longest).min(max_distance) {
            return EdgeSnap::Far;
        }
        // On the side opposite corner k: halfedge k + 1, from corner k + 1 to k + 2.
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        let h = 3 * face as usize + i;
        let t_halfedge = bary[j] / (bary[i] + bary[j]);
        if !(t_halfedge > 0.0 && t_halfedge < 1.0) {
            return EdgeSnap::Corner;
        }
        let edge = mesh.edge(h);
        let forward = mesh.from_vertex(h) == mesh.edge_vertices(edge)[0];
        // The direction of the edge in the frame of the face (its halfedge 0).
        let d = [v[j][0] - v[i][0], v[j][1] - v[i][1]];
        let along = d[1].atan2(d[0]) + if forward { 0.0 } else { PI };
        EdgeSnap::Onto {
            location: SurfacePoint::Edge {
                edge: edge as u32,
                t: if forward {
                    t_halfedge
                } else {
                    1.0 - t_halfedge
                },
            },
            back: wrap(back - along),
            coordinate: bary[k],
        }
    }

    /// The corner angle at the tail of halfedge `h`, which exists: the caller has checked the
    /// face.
    pub(crate) fn angle_at(&self, h: usize) -> f64 {
        let l = self.face_lengths(h / 3);
        let k = h % 3;
        corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3]).unwrap_or(0.0)
    }

    /// The signpost of the next halfedge counter-clockwise after `h` around its tail, from
    /// the signpost of `h` and the corner angle between them.
    pub(crate) fn signpost_after(&self, h: usize) -> f64 {
        let v = self.from_vertex(h);
        signpost_update(
            self.signposts[h],
            self.angle_at(h),
            self.angle_sums[v],
            self.is_boundary_vertex(v),
        )
    }

    /// Inserts a vertex at the point `bary` of face `f` (barycentric coordinates, all
    /// positive) and joins it to the three corners (paper Alg. 7). The face keeps its id and
    /// becomes the triangle on its first side; two faces, three edges and one vertex are
    /// appended.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`] for a face that does not exist or coordinates that are
    /// not positive numbers summing to one, and [`IntrinsicError::InsertionRefused`] when a
    /// new triangle would be too thin (the point is next to an edge of the face: split that
    /// edge instead) or the place of the point on the input mesh cannot be found. The
    /// triangulation is unchanged in every case.
    pub fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntrinsicError> {
        let at = SurfacePoint::Face {
            face: f as u32,
            bary,
        };
        if let Some(reason) = point_defect(&at, self.n_vertices(), self.n_edges(), self.n_faces()) {
            return Err(IntrinsicError::InvalidPoint { point: at, reason });
        }
        let refused = |reason| IntrinsicError::InsertionRefused { at, reason };
        let l = self.face_lengths(f);
        let layout = triangle_layout(l).ok_or(IntrinsicError::DegenerateFace {
            surface: Self::NAME,
            face: f as u32,
            lengths: l,
        })?;
        let d = corner_distances(l, bary).ok_or(IntrinsicError::DegenerateFace {
            surface: Self::NAME,
            face: f as u32,
            lengths: l,
        })?;
        // The three triangles: on side k, with the new edges d[k] and d[k + 1]. Twice the area
        // of the triangle on side k is bary[k + 2] times twice the area of the face.
        let twice_area = layout[1][0] * layout[2][1];
        let mut shape = f64::INFINITY;
        let mut thinnest = 0;
        for k in 0..3 {
            let (a, b) = (d[k], d[(k + 1) % 3]);
            let longest = l[k].max(a).max(b);
            let s = bary[(k + 2) % 3] * twice_area / (a.max(b) * longest);
            if !triangle_is_valid(l[k], b, a) {
                shape = 0.0;
                thinnest = k;
            } else if s < shape {
                shape = s;
                thinnest = k;
            }
        }
        if shape < FLIP_SHAPE_LIMIT || shape.is_nan() {
            return Err(refused(InsertRefusal::TooThin {
                shape,
                nearest_edge: Some(self.complex.halfedge_edge[3 * f + thinnest]),
            }));
        }
        // The place on the input mesh, traced from the nearest corner.
        let x = (0..3).min_by(|&a, &b| d[a].total_cmp(&d[b])).unwrap_or(0);
        let hx = 3 * f + x;
        let vx = self.complex.faces[f][x] as usize;
        let (u, w) = (
            [
                layout[(x + 1) % 3][0] - layout[x][0],
                layout[(x + 1) % 3][1] - layout[x][1],
            ],
            {
                let p = [
                    bary[1] * layout[1][0] + bary[2] * layout[2][0],
                    bary[1] * layout[1][1] + bary[2] * layout[2][1],
                ];
                [p[0] - layout[x][0], p[1] - layout[x][1]]
            },
        );
        let inside = (u[0] * w[1] - u[1] * w[0]).atan2(u[0] * w[0] + u[1] * w[1]);
        let phi = signpost_update(
            self.signposts[hx],
            inside.max(0.0),
            self.angle_sums[vx],
            self.is_boundary_vertex(vx),
        );
        let placed = self.place_by_trace(vx, phi, d[x]).map_err(refused)?;

        // Connectivity (see `FaceSplit`). Before: f = (a, b, c). After: f = (a, b, p),
        // f1 = (b, c, p), f2 = (c, a, p); b → c and c → a move to f1 and f2 with their
        // signposts.
        let split = self.complex.split_face_unchecked(f);
        split.move_halfedge_data(&mut self.signposts, 0.0);
        let p = split.vertex;
        let [f1, f2] = split.new_faces.map(|g| g as usize);
        let [ea, eb, ec] = split.new_edges;
        // The new edges are appended in this order, from corners a, b, c.
        self.edge_lengths.extend(d);
        self.edge_on_input.extend([NONE; 3]);
        self.angle_sums.push(TAU);
        self.locations.push(placed.location);

        // Signposts of the halfedges that arrive: each follows the halfedge of its corner.
        self.signposts[3 * f2 + 1] = self.signpost_after(3 * f);
        self.signposts[3 * f + 1] = self.signpost_after(3 * f1);
        self.signposts[3 * f1 + 1] = self.signpost_after(3 * f2);
        // Signposts of the halfedges that leave, counter-clockwise: to a, to b, to c.
        let out = [3 * f + 2, 3 * f1 + 2, 3 * f2 + 2];
        self.signposts[out[x]] = placed.back;
        for step in 0..2 {
            let from = out[(x + step) % 3];
            self.signposts[out[(x + step + 1) % 3]] = self.signpost_after(from);
        }
        Ok(Insertion {
            vertex: p,
            location: placed.location,
            new_edges: vec![ea, eb, ec],
            new_faces: vec![f1 as u32, f2 as u32],
            snapped: placed.snapped,
        })
    }

    /// Splits edge `e` at the parameter `t` from the tail of its halfedge (`0 < t < 1`): a
    /// vertex is inserted there and joined to the one or two opposite vertices (paper Alg. 8).
    /// The half from the tail keeps the id `e`; the faces of the edge keep their ids and
    /// become the triangles on that half.
    ///
    /// A vertex on an edge that lies along an input edge ([`edge_on_input`]), which every
    /// boundary edge does, is placed on that input edge by interpolation, and both halves
    /// still lie along it. Otherwise its place is traced from the nearer end.
    ///
    /// [`edge_on_input`]: IntrinsicTriangulation::edge_on_input
    ///
    /// # Errors
    /// [`IntrinsicError::IndexOutOfRange`], [`IntrinsicError::InvalidPoint`] for a parameter
    /// that is not strictly between 0 and 1, and [`IntrinsicError::InsertionRefused`] as for
    /// [`insert_vertex`](Self::insert_vertex). The triangulation is unchanged in every case.
    pub fn split_edge(&mut self, e: usize, t: f64) -> Result<Insertion, IntrinsicError> {
        self.check_edge(e)?;
        let at = SurfacePoint::Edge { edge: e as u32, t };
        if !(t > 0.0 && t < 1.0) {
            return Err(IntrinsicError::InvalidPoint {
                point: at,
                reason: "the edge parameter is not strictly between 0 and 1",
            });
        }
        let refused = |reason| IntrinsicError::InsertionRefused { at, reason };
        let h0 = self.complex.edge_halfedge[e] as usize;
        let (f0, k0) = (h0 / 3, h0 % 3);
        let h1 = self.complex.twin[h0];
        if h1 != NONE && h1 as usize / 3 == f0 {
            return Err(refused(InsertRefusal::SameFaceOnBothSides));
        }
        let length = self.edge_lengths[e];
        let (la, lb) = (t * length, (1.0 - t) * length);
        let degenerate = |tri: &Self, f: usize| IntrinsicError::DegenerateFace {
            surface: Self::NAME,
            face: f as u32,
            lengths: tri.face_lengths(f),
        };
        // The distance to the opposite vertex of a side: the edge is side `k` of face `f`,
        // the point at parameter `s` from the tail of that side.
        let across = |tri: &Self, f: usize, k: usize, s: f64| {
            let l = tri.face_lengths(f);
            let l = [l[k], l[(k + 1) % 3], l[(k + 2) % 3]];
            corner_distances(l, [1.0 - s, s, 0.0])
                .map(|d| (d[2], l))
                .ok_or_else(|| degenerate(tri, f))
        };
        // Shape of the two triangles a side is cut into: (tail, p, apex) and (p, head, apex).
        let shapes = |tail: f64, head: f64, to_apex: f64, l: [f64; 3]| {
            let twice_area = {
                let v = triangle_layout(l).unwrap_or([[0.0; 2]; 3]);
                v[1][0] * v[2][1]
            };
            let (s, r) = (tail / l[0], head / l[0]);
            let first = if triangle_is_valid(tail, to_apex, l[2]) {
                s * twice_area / (tail.max(to_apex) * tail.max(to_apex).max(l[2]))
            } else {
                0.0
            };
            let second = if triangle_is_valid(head, l[1], to_apex) {
                r * twice_area / (head.max(to_apex) * head.max(to_apex).max(l[1]))
            } else {
                0.0
            };
            first.min(second)
        };
        let (dc, l0) = across(self, f0, k0, t)?;
        let mut shape = shapes(la, lb, dc, l0);
        let mut dd = 0.0;
        if h1 != NONE {
            let h1 = h1 as usize;
            let (d, l1) = across(self, h1 / 3, h1 % 3, 1.0 - t)?;
            dd = d;
            shape = shape.min(shapes(lb, la, d, l1));
        }
        if shape < FLIP_SHAPE_LIMIT || shape.is_nan() {
            return Err(refused(InsertRefusal::TooThin {
                shape,
                nearest_edge: None,
            }));
        }
        let a = self.complex.faces[f0][k0] as usize;
        let b = self.complex.faces[f0][(k0 + 1) % 3] as usize;
        // The place of the new vertex on the input mesh, and the signpost of one halfedge
        // that leaves it: `Ok(direction to b)` or `Err(direction to a)`.
        let on_input = self.edge_on_input[e];
        let along = if on_input == NONE {
            None
        } else {
            match (
                self.parameter_on_input_edge(a, on_input as usize),
                self.parameter_on_input_edge(b, on_input as usize),
            ) {
                (Some(ta), Some(tb)) if ta != tb => Some((ta, tb)),
                _ => None,
            }
        };
        let (location, first, snapped) = match along {
            Some((ta, tb)) => {
                let location = SurfacePoint::Edge {
                    edge: on_input,
                    t: (1.0 - t) * ta + t * tb,
                };
                // On the boundary the frame starts at the halfedge that leaves towards b.
                let to_b = if h1 == NONE || tb > ta { 0.0 } else { PI };
                (location, Ok(to_b), None)
            }
            None if h1 == NONE => {
                return Err(refused(InsertRefusal::NotLocated {
                    detail: format!("boundary edge {e} is not known to lie along an input edge"),
                }));
            }
            None => {
                let placed = if t <= 0.5 {
                    self.place_by_trace(a, self.signposts[h0], la)
                } else {
                    self.place_by_trace(b, self.signposts[h1 as usize], lb)
                }
                .map_err(refused)?;
                let first = if t <= 0.5 {
                    Err(placed.back)
                } else {
                    Ok(placed.back)
                };
                (placed.location, first, placed.snapped)
            }
        };

        // Connectivity (see `EdgeSplit`). Before: f0 = (a, b, c) with a → b in slot k0, and
        // f1 = (b, a, d). After: f0 = (a, p, c), f2 = (p, b, c), f1 = (p, a, d),
        // f3 = (b, p, d); b → c, d → b and b → a (as b → p) move with their signposts.
        let split = self.complex.split_edge_unchecked(e);
        split.move_halfedge_data(&mut self.signposts, 0.0);
        let p = split.vertex;
        let e2 = split.halves[1] as usize;
        let [ec, ed] = split.to_opposite;
        let [f2, f3] = split.new_faces;
        let f2 = f2 as usize;
        self.edge_lengths.resize(self.complex.n_edges(), 0.0);
        self.edge_on_input.resize(self.complex.n_edges(), NONE);
        self.edge_lengths[e2] = lb;
        self.edge_on_input[e2] = on_input;
        if ed != NONE {
            self.edge_lengths[ed as usize] = dd;
        }
        self.edge_lengths[ec as usize] = dc;
        self.edge_lengths[e] = la;
        self.angle_sums.push(if h1 == NONE { PI } else { TAU });
        self.locations.push(location);
        let out: Vec<usize> = split.outgoing().iter().map(|&h| h as usize).collect();
        let new_edges = split.new_edges();
        let new_faces: Vec<u32> = split.new_faces.into_iter().filter(|&g| g != NONE).collect();

        // Signposts of the halfedges that arrive from the opposite vertices.
        self.signposts[3 * f2 + 2] = self.signpost_after(3 * f0 + (k0 + 2) % 3);
        if h1 != NONE {
            let h1 = h1 as usize;
            let f3 = f3 as usize;
            self.signposts[3 * (h1 / 3) + (h1 % 3 + 2) % 3] = self.signpost_after(3 * f3 + 2);
        }
        // Signposts of the halfedges that leave, counter-clockwise: to b, to c, to a, to d.
        let start = match first {
            Ok(to_b) => {
                self.signposts[out[0]] = to_b;
                0
            }
            Err(to_a) => {
                self.signposts[out[2]] = to_a;
                2
            }
        };
        let n = out.len();
        for step in 0..n - 1 {
            let from = out[(start + step) % n];
            self.signposts[out[(start + step + 1) % n]] = self.signpost_after(from);
        }
        Ok(Insertion {
            vertex: p,
            location,
            new_edges,
            new_faces,
            snapped,
        })
    }
}
