//! Tracing over a triangulation that stores no signposts (SPEC N50; IC-SPEC I9).
//!
//! The tracing kernel ([`crate::trace`]) needs, at a vertex, the direction in which every
//! halfedge leaves it (its signpost) and the angle sum. A representation that stores only the
//! connectivity and the edge lengths (the integer coordinates of Gillespie, Sharp and Crane
//! 2021) has both implicitly: lay the corners around the vertex out from the lengths.
//! [`AngleFrames`] does that on demand and hands the result to the same kernel that traces over
//! [`IntrinsicTriangulation`](crate::IntrinsicTriangulation) and the input mesh.
//!
//! # The frame of a vertex
//!
//! Directions at a vertex `v` are normalized angles as in the signpost representation (paper
//! Sec. 3.1, SPEC N2): the corner angles around `v` are scaled by `2π / Θ` at an interior vertex
//! and by `π / Θ` at a boundary vertex, with `Θ` here the sum of the corner angles computed
//! from the lengths. The direction `0` is that of [`Complex::outgoing_halfedge`]: the boundary
//! halfedge leaving a boundary vertex (so the frame of a boundary vertex is the stored one), an
//! arbitrary but fixed halfedge at an interior vertex. The signpost of the `k`-th outgoing
//! halfedge counter-clockwise is `s · (θ_0 + … + θ_{k−1})` (`it-angle-normalization`; the
//! running sum is scaled once, as for the input mesh, instead of being updated corner by
//! corner, so it stays in `[0, 2π)` without the reduction of SPEC N18).
//!
//! A frame belongs to the current connectivity: a flip of an edge at an interior vertex can
//! change which halfedge is [`Complex::outgoing_halfedge`], and with it the frame. Directions
//! are to be used right after they are computed (a move of a vertex traces at once, IC-SPEC
//! §4.7), not stored across operations.
//!
//! # Agreement with stored signposts
//!
//! On an [`IntrinsicTriangulation`](crate::IntrinsicTriangulation) the computed signpost of a
//! halfedge equals the stored one minus the stored signpost of the vertex's outgoing halfedge
//! (reduced to `[0, 2π)`), up to the rounding of the stored updates; `tests/angle_frames.rs`
//! measures the difference after random flips and insertions and compares traces.

use std::cell::{Ref, RefCell};
use std::f64::consts::TAU;

use ddg_core::geometry::corner_angle;

use crate::complex::Complex;
use crate::error::IntrinsicError;
use crate::formulas::angle_scale;
use crate::interface::{Triangulation, face_lengths};
use crate::point::SurfacePoint;
use crate::surface::Surface;
use crate::trace::{Trace, TraceOptions, trace_query};

/// A view of a triangulation (connectivity and edge lengths) whose signposts and angle sums
/// are computed from the corner angles when the tracing kernel asks for them. See the
/// [module documentation](self) for the frames.
///
/// The fan of the vertex asked last is kept, so the kernel's repeated questions about one
/// vertex cost one layout of its corners.
#[derive(Debug)]
pub struct AngleFrames<'a> {
    complex: &'a Complex,
    lengths: &'a [f64],
    fan: RefCell<Fan>,
}

/// The computed directions around one vertex.
#[derive(Debug, Default)]
struct Fan {
    /// The vertex, or `usize::MAX` before the first question.
    vertex: usize,
    /// The sum of its corner angles.
    theta: f64,
    /// `(halfedge, signpost)` for its outgoing halfedges, sorted by halfedge.
    entries: Vec<(u32, f64)>,
}

impl<'a> AngleFrames<'a> {
    /// The view of a connectivity with edge lengths.
    ///
    /// # Errors
    /// [`IntrinsicError::Inconsistent`] if there is not one length per edge.
    pub fn new(
        complex: &'a Complex,
        lengths: &'a [f64],
    ) -> Result<AngleFrames<'a>, IntrinsicError> {
        if lengths.len() != complex.n_edges() {
            return Err(IntrinsicError::Inconsistent {
                what: format!(
                    "{} edge lengths for {} edges",
                    lengths.len(),
                    complex.n_edges()
                ),
            });
        }
        Ok(AngleFrames {
            complex,
            lengths,
            fan: RefCell::new(Fan {
                vertex: usize::MAX,
                ..Fan::default()
            }),
        })
    }

    /// The view of any representation.
    pub fn of<T: Triangulation>(triangulation: &'a T) -> AngleFrames<'a> {
        AngleFrames {
            complex: triangulation.complex(),
            lengths: triangulation.edge_lengths(),
            fan: RefCell::new(Fan {
                vertex: usize::MAX,
                ..Fan::default()
            }),
        }
    }

    /// The fan of vertex `v`, computed if it is not the one kept.
    fn fan(&self, v: usize) -> Ref<'_, Fan> {
        if self.fan.borrow().vertex != v {
            let mut fan = self.fan.borrow_mut();
            let c = self.complex;
            fan.vertex = v;
            fan.entries.clear();
            let mut theta = 0.0;
            // First the corner angles, then their running sums scaled.
            for h in c.outgoing_halfedges(v) {
                let l = face_lengths(c, self.lengths, h / 3);
                let k = h % 3;
                let angle = corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3]).unwrap_or(0.0);
                fan.entries.push((h as u32, angle));
                theta += angle;
            }
            let scale = angle_scale(theta, c.is_boundary_vertex(v));
            let mut running = 0.0;
            for entry in fan.entries.iter_mut() {
                let angle = entry.1;
                let phi = running * scale;
                entry.1 = if phi >= TAU { phi - TAU } else { phi };
                running += angle;
            }
            fan.theta = theta;
            fan.entries.sort_unstable_by_key(|entry| entry.0);
        }
        self.fan.borrow()
    }

    /// The direction in which halfedge `h` leaves its tail, in the frame of the tail (see the
    /// [module documentation](self)).
    ///
    /// # Panics
    /// If `h` is out of range.
    pub fn signpost(&self, h: usize) -> f64 {
        let v = self.complex.from_vertex(h);
        let fan = self.fan(v);
        match fan
            .entries
            .binary_search_by_key(&(h as u32), |entry| entry.0)
        {
            Ok(i) => fan.entries[i].1,
            // Every halfedge leaving `v` is in its fan in a valid complex.
            Err(_) => f64::NAN,
        }
    }

    /// The sum of the corner angles at vertex `v`, from the lengths.
    ///
    /// # Panics
    /// If `v` is out of range.
    pub fn angle_sum(&self, v: usize) -> f64 {
        self.fan(v).theta
    }

    /// Follows the straight path that leaves `start` (a point of the triangulation) in
    /// `direction` for `length` over the triangulation, as
    /// [`IntrinsicTriangulation::trace`](crate::IntrinsicTriangulation::trace) does with
    /// stored signposts: from a vertex, `direction` is in the frame of this view; on an edge or
    /// in a face, in the frame of the edge or face (SPEC N7, [`SurfacePoint`]).
    ///
    /// # Errors
    /// As [`InputSurface::trace`](crate::InputSurface::trace).
    pub fn trace(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
        options: &TraceOptions,
    ) -> Result<Trace, IntrinsicError> {
        trace_query(self, start, direction, length, options)
    }
}

impl Surface for AngleFrames<'_> {
    const NAME: &'static str = "intrinsic triangulation";

    #[inline]
    fn n_vertices(&self) -> usize {
        self.complex.n_vertices()
    }
    #[inline]
    fn n_edges(&self) -> usize {
        self.complex.n_edges()
    }
    #[inline]
    fn n_faces(&self) -> usize {
        self.complex.n_faces()
    }
    #[inline]
    fn twin(&self, h: usize) -> u32 {
        self.complex.twin[h]
    }
    #[inline]
    fn tail(&self, h: usize) -> u32 {
        self.complex.faces[h / 3][h % 3]
    }
    #[inline]
    fn edge(&self, h: usize) -> u32 {
        self.complex.halfedge_edge[h]
    }
    #[inline]
    fn length(&self, e: usize) -> f64 {
        self.lengths[e]
    }
    #[inline]
    fn outgoing(&self, v: usize) -> u32 {
        self.complex.vertex_halfedge[v]
    }
    fn signpost(&self, h: usize) -> f64 {
        AngleFrames::signpost(self, h)
    }
    fn angle_sum(&self, v: usize) -> f64 {
        AngleFrames::angle_sum(self, v)
    }
    #[inline]
    fn is_forward(&self, h: usize) -> bool {
        self.complex.is_forward(h)
    }
    #[inline]
    fn edge_halfedge(&self, e: usize) -> u32 {
        self.complex.edge_halfedge[e]
    }
}
