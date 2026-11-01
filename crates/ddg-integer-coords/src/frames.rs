//! Input-aligned tangent frames on integer coordinates (SPEC §4.12, I26).
//!
//! The signposts of `ddg_intrinsic` store, per halfedge, its direction at its tail in a frame
//! that is the input mesh's at an input vertex and that of the input face or edge at an inserted
//! vertex (IT-SPEC N2, N7, N9, N10). Integer coordinates store no angles, but the same
//! directions follow from what they do store: the roundabout of a halfedge at an input vertex
//! names the input wedge it leaves into, and the straight segment of its edge over the input mesh
//! (`ExtractIntrinsicEdge`, SPEC §3.1, whose input faces the integers decide) gives the angle
//! inside that wedge, or inside the input face or on the input edge that an inserted vertex lies
//! in. [`IntegerTriangulation::frame_directions`] computes them for every halfedge in one pass;
//! on the same triangulation they equal the signposts' `signpost(h)` to rounding.
//!
//! The frames, as the signposts define them (normalized, so that a full turn is `2π` inside and
//! `π` on the boundary):
//!
//! * **an input vertex `a`**: angle 0 along the input halfedge of signpost 0
//!   (`InputSurface::signpost`), the other input halfedges at the sum of the corner angles before
//!   them times `s_a = 2π/Θ_a` (`π/Θ_a` on the boundary). A halfedge `h` with roundabout `r`
//!   leaving along an input edge (`n = −1`) has the input signpost of the input halfedge of local
//!   index `r` (`π` for the arriving boundary edge); one leaving into the wedge between local
//!   indices `r − 1` and `r` has `φ_{r−1} + s_a α`, `α` the angle from the input halfedge `r − 1`
//!   to the segment, in the layout of that wedge's input face;
//! * **an inserted vertex in input face `g`**: angle 0 along the halfedge `3g`, unscaled:
//!   `φ = arg(d)` measured from that side in the layout of `g`;
//! * **an inserted vertex on interior input edge `E`**: angle 0 along `E` towards its
//!   higher-numbered vertex (`edge_vertices(E)[1]`); in the input face of the halfedge `x` of `E`
//!   the segment leaves into, `φ = β_x + α`, `β_x = 0` if `x` runs from `edge_vertices(E)[0]`,
//!   else `π`, `α ∈ [0, π]` the angle from `x`; a halfedge along `E` has `0` towards the larger
//!   parameter, `π` towards the smaller;
//! * **an inserted vertex on a boundary input edge**: angle 0 along the boundary halfedge that
//!   leaves it (the input boundary halfedge, the surface on its left), `π` along the arriving
//!   one: `φ = α` from the input boundary halfedge.
//!
//! Angles are reduced to `[0, 2π)` inside and kept in `[0, π]` on the boundary, as the signposts
//! store them. The angle `α` inside a corner is clamped to `[0, corner]` (counted, SPEC I5: zero
//! with exact arithmetic), as the signposts' tracing clamps its back direction.
//!
//! **Conditioning.** `α` comes from the direction of the whole straight segment of the edge in
//! its unfolded strip, so its error is the rounding of the strip's layout divided by the edge's
//! length, not by the length of its first piece: an edge that leaves an inserted vertex next to
//! an input edge is as well conditioned as any other.
//!
//! **Cost.** One [`extract_all_edges`](IntegerTriangulation::extract_all_edges) (which labels
//! every crossing with its input edge), then per edge with `n ≥ 0` the layout of its strip of
//! input faces: `O(|E⁰| + |E¹| + Σ n⁺)` in all, without the walks of `ExtractIntrinsicEdge`.

use std::f64::consts::{PI, TAU};

use ddg_intrinsic::SurfacePoint;
use ddg_intrinsic::formulas::angle_scale;

use crate::coords::IntegerTriangulation;
use crate::error::{IntegerError, Place, ValidityItem};
use crate::geometry::{ChordEnd, ChordEndPlace};
use crate::trace::PathPiece;

/// What [`IntegerTriangulation::frame_directions`] did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FrameReport {
    /// Halfedges whose edge runs along an input edge (`n = −1`): the direction is read off the
    /// integers alone.
    pub n_along: usize,
    /// Halfedges whose direction comes from the straight segment of their edge.
    pub n_from_segment: usize,
    /// Directions inside a corner of an input vertex whose angle `α` fell outside the corner by
    /// rounding and was clamped (SPEC I5).
    pub n_clamped: usize,
    /// The largest clamp, in radians.
    pub max_clamp: f64,
    /// The largest relative difference between a segment's length in its unfolded strip and the
    /// edge's length (as `StripReport::max_relative_mismatch`).
    pub max_relative_mismatch: f64,
}

fn invalid(h: usize, detail: String) -> IntegerError {
    IntegerError::invalid(ValidityItem::Extraction, Place::Halfedge(h as u32), detail)
}

/// `angle` reduced to `[0, 2π)`.
fn wrap(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    if a >= TAU { 0.0 } else { a }
}

/// The angle from `u` to `d`, counter-clockwise, in `(−π, π]`.
fn angle_from(u: [f64; 2], d: [f64; 2]) -> f64 {
    (u[0] * d[1] - u[1] * d[0]).atan2(u[0] * d[0] + u[1] * d[1])
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

impl IntegerTriangulation {
    /// The direction of every halfedge at its tail in the input-aligned frame of the tail,
    /// normalized as the signposts store it (`IntrinsicTriangulation::signpost`: `[0, 2π)` at an
    /// interior vertex, `[0, π]` at a boundary vertex; the [module documentation](self) gives the
    /// frames). On the same triangulation (connectivity, lengths, places of the inserted vertices)
    /// these are the signposts to rounding.
    ///
    /// [formula: ic-input-aligned-frame]
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] naming the halfedge, edge or crossing when the
    /// integers are not those of a realization (what `certify` rejects), and the errors of
    /// [`extract_all_edges`](Self::extract_all_edges) and of the strip layout (a degenerate
    /// input face).
    pub fn frame_directions(&self) -> Result<(Vec<f64>, FrameReport), IntegerError> {
        let complex = &self.complex;
        let ne = complex.n_edges();
        // The input edge of every crossing, per edge from the tail of its edge_halfedge.
        let mut offsets = Vec::with_capacity(ne + 1);
        let mut total = 0usize;
        for e in 0..ne {
            offsets.push(total);
            total += self.normal[e].max(0) as usize;
        }
        offsets.push(total);
        let mut labels = vec![u32::MAX; total];
        for path in self.extract_all_edges()? {
            for piece in &path.pieces {
                let PathPiece::Crossings(segment) = piece else {
                    continue;
                };
                for z in &segment.crossings {
                    let h = z.halfedge as usize;
                    let e = complex.edge(h);
                    let n = self.normal[e];
                    if z.index < 0 || z.index >= n {
                        return Err(invalid(
                            h,
                            format!("crossing {} of an edge with n = {n}", z.index),
                        ));
                    }
                    let p = if complex.is_forward(h) {
                        z.index
                    } else {
                        n - 1 - z.index
                    } as usize;
                    let slot = &mut labels[offsets[e] + p];
                    if *slot != u32::MAX {
                        return Err(invalid(
                            h,
                            format!(
                                "crossing {p} of edge {e} lies on input edges {} and {}",
                                *slot, path.input_edge
                            ),
                        ));
                    }
                    *slot = path.input_edge;
                }
            }
        }
        if let Some(i) = labels.iter().position(|&l| l == u32::MAX) {
            let e = offsets.partition_point(|&o| o <= i) - 1;
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Edge(e as u32),
                format!("crossing {} lies on no input edge", i - offsets[e]),
            ));
        }
        let mut out = vec![0.0; complex.n_halfedges()];
        let mut report = FrameReport::default();
        for e in 0..ne {
            let h = complex.edge_halfedge(e);
            let twin = complex.twin(h);
            if self.normal[e] < 0 {
                out[h] = self.along_direction(h)?;
                report.n_along += 1;
                if let Some(t) = twin {
                    out[t] = self.along_direction(t)?;
                    report.n_along += 1;
                }
                continue;
            }
            let (chord, _) = self.input_chord(e, Some(&labels[offsets[e]..offsets[e + 1]]))?;
            report.max_relative_mismatch =
                report.max_relative_mismatch.max(chord.relative_mismatch);
            out[h] = self.segment_direction(h, &chord.tail, &chord.head, &mut report)?;
            report.n_from_segment += 1;
            let Some(t) = twin else {
                return Err(invalid(
                    h,
                    format!("boundary edge {e} has n = {}", self.normal[e]),
                ));
            };
            out[t] = self.segment_direction(t, &chord.head, &chord.tail, &mut report)?;
            report.n_from_segment += 1;
        }
        Ok((out, report))
    }

    /// The direction of halfedge `h` (`n = −1`) at its tail: along an input edge.
    fn along_direction(&self, h: usize) -> Result<f64, IntegerError> {
        let complex = &self.complex;
        let mesh = self.input.mesh();
        let a = complex.from_vertex(h);
        if self.is_input_vertex(a) {
            let r = self.roundabout(h).ok_or_else(|| {
                invalid(h, format!("leaves input vertex {a} without a roundabout"))
            })?;
            let deg = self.local.degree(a);
            if r >= deg {
                return Err(invalid(
                    h,
                    format!("roundabout {r} at vertex {a} is not below deg₀ = {deg}"),
                ));
            }
            if self.local.is_boundary(a) && r == deg - 1 {
                // The boundary edge that arrives at `a`: the constant `π` (IT-SPEC N2).
                return Ok(PI);
            }
            let edge = self.local.edge(a, r);
            let x = mesh.edge_halfedge(edge);
            let out = if mesh.from_vertex(x) == a {
                Some(x)
            } else {
                mesh.twin(x).filter(|&t| mesh.from_vertex(t) == a)
            }
            .ok_or_else(|| {
                invalid(
                    h,
                    format!("input edge {edge} has no halfedge leaving vertex {a}"),
                )
            })?;
            return Ok(self.input.signpost(out));
        }
        let SurfacePoint::Edge { edge, t } = self.locations[a] else {
            return Err(invalid(
                h,
                format!(
                    "runs along an input edge, but its tail {a} lies at {}",
                    self.locations[a]
                ),
            ));
        };
        if complex.twin(h).is_none() {
            // The boundary halfedge that leaves an inserted boundary vertex.
            return Ok(0.0);
        }
        let b = complex.to_vertex(h);
        let tb = self.parameter_on_input_edge(b, edge as usize)?;
        Ok(if tb > t { 0.0 } else { PI })
    }

    /// The direction of halfedge `h` at its tail, from the straight segment of its edge: `from`
    /// is the end at the tail, `to` the other end.
    fn segment_direction(
        &self,
        h: usize,
        from: &ChordEnd,
        to: &ChordEnd,
        report: &mut FrameReport,
    ) -> Result<f64, IntegerError> {
        let a = self.complex.from_vertex(h);
        let mesh = self.input.mesh();
        let c = &from.corners;
        let d = sub(to.point, from.point);
        let boundary = self.complex.is_boundary_vertex(a);
        match from.place {
            ChordEndPlace::Corner { out } => {
                if !self.is_input_vertex(a) || mesh.from_vertex(out) != a {
                    return Err(invalid(
                        h,
                        format!(
                            "its segment starts at a corner of input halfedge {out}, not at vertex {a}"
                        ),
                    ));
                }
                let s = out % 3;
                let corner = self.input.geometry().corner_angles()[out / 3][s];
                let raw = angle_from(sub(c[(s + 1) % 3], c[s]), d);
                let alpha = raw.clamp(0.0, corner);
                let clamp = (raw - alpha).abs();
                if clamp > 0.0 {
                    report.n_clamped += 1;
                    report.max_clamp = report.max_clamp.max(clamp);
                }
                let theta = self.input.vertex_angle_sum(a);
                let phi = self.input.signpost(out) + angle_scale(theta, boundary) * alpha;
                Ok(if boundary { phi.min(PI) } else { wrap(phi) })
            }
            ChordEndPlace::OnEdge { side } => {
                let s = side % 3;
                let alpha = angle_from(sub(c[(s + 1) % 3], c[s]), d);
                let edge = mesh.edge(side);
                if boundary {
                    // The input boundary halfedge is the side itself: the frame's direction 0.
                    Ok(alpha.clamp(0.0, PI))
                } else {
                    let base = if mesh.from_vertex(side) == mesh.edge_vertices(edge)[0] {
                        0.0
                    } else {
                        PI
                    };
                    Ok(wrap(base + alpha))
                }
            }
            ChordEndPlace::InFace => Ok(wrap(angle_from(sub(c[1], c[0]), d))),
        }
    }
}
