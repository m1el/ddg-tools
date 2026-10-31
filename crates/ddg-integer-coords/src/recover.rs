//! `RecoverBarycentric` (paper App. A, Eq. 15–16; SPEC §4.4, I16) and the point query
//! `T¹ → T⁰` (SPEC §4.9).
//!
//! A region `R` of a face of `T¹` (the face cut by the arcs of the input edges that cross it,
//! [`crate::insert`]) lies inside one input face `abc`: no input edge crosses its interior. Its
//! corners (3 to 6: corners of the face and crossings on its sides) are known on both sides, as
//! barycentric coordinates `v^(m)` on the face of `T¹` and `u^(m)` on `abc` (a crossing is a
//! point of its input edge, a corner of the face is the place `q⁰` of its vertex). A point `v`
//! of the region is the affine combination `Σ ξ_m v^(m)` (Eq. 15; the minimum-norm `ξ`, by QR of
//! the `ρ × 3` matrix of the corners, SPEC I16) and its place is `Σ ξ_m u^(m)` (Eq. 16). With exact
//! corners this is the point exactly, for any solution `ξ` (Lean `recover_eq`); every corner is
//! used (the reference picks three or four, map §5 item 2).
//!
//! **The input face** of the region is read off the integers: every region meets the boundary of
//! the face in a gap between two consecutive crossings of a side (Lean `FaceSplit.lean`, regions
//! as boundary gaps), and `ExtractIntrinsicEdge` of that side says in which input face the side
//! runs there (the face its path leaves through the next crossing, or enters through the
//! previous one, or its first face for a side without crossings). A face whose three sides run
//! along input edges is an input face; it is found from the roundabout of a side at an input
//! vertex. Every corner of the region must be a point of that input face's closure, or the
//! recovery fails ([`IntegerError::RecoverFailed`]).
//!
//! **Floating point decides only the position.** The coordinates on `abc` are reported as
//! computed ([`RecoveredPoint`], with the smallest one before anything is done to it). For an
//! inserted vertex the place must be a point strictly inside `abc` (SPEC §2.2): a point whose
//! coordinates are within [`RECOVER_ZERO`] of the range `(0, 1)` is moved, by the smallest
//! amount, towards the centroid of the region's corners on `abc` (which lies inside the region)
//! until every coordinate is at least `RECOVER_ZERO` times that of the centroid, and the amount is
//! reported ([`RecoverReport::moved`]); farther outside is [`IntegerError::RecoverFailed`]. The
//! point query instead returns a point within rounding of a side (or corner) of `abc` as a point
//! of that input edge (or that input vertex), as the signposts' queries do.

use ddg_intrinsic::{EdgeTraceOptions, SurfacePoint};

use crate::coords::IntegerTriangulation;
use crate::error::{IntegerError, Place, ValidityItem};
use crate::insert::{FaceArcs, Region, RegionPoint};

/// Barycentric coordinates on the input face that are below `0` or above `1` by at most this
/// are rounding (SPEC I16); farther outside is [`IntegerError::RecoverFailed`].
pub const RECOVER_ZERO: f64 = 1e-12;

/// App. A, Eq. 15–16: the barycentric coordinates on the input face of the point with
/// coordinates `v` on the face of the triangulation, from the corners of its region with
/// coordinates `on_triangulation[m]` on the face of `T¹` and `on_input[m]` on the input face.
/// `ξ` is the minimum-norm solution of the `3 × ρ` system `Σ ξ_m on_triangulation[m] = v`
/// (Eq. 15), computed by a QR factorization of its transpose (modified Gram–Schmidt with one
/// reorthogonalization; `ρ ≤ 6`), and the result is `Σ ξ_m on_input[m]` (Eq. 16). Returns the
/// coordinates and the residual `max_i |Σ_m ξ_m v^(m)_i − v_i|`; `None` for fewer than three
/// corners, slices of different lengths, or corners that do not span the plane (a zero pivot).
///
/// [formula: ic-recover-barycentric]
pub fn recover_barycentric(
    on_triangulation: &[[f64; 3]],
    on_input: &[[f64; 3]],
    v: [f64; 3],
) -> Option<([f64; 3], f64)> {
    let rho = on_triangulation.len();
    if rho < 3 || on_input.len() != rho {
        return None;
    }
    // The columns of Bᵀ (ρ × 3): column i is the i-th coordinate of every corner.
    let mut q: Vec<[f64; 3]> = on_triangulation.to_vec();
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for _pass in 0..2 {
            for j in 0..i {
                let d: f64 = q.iter().map(|row| row[j] * row[i]).sum();
                r[j][i] += d;
                for row in q.iter_mut() {
                    row[i] -= d * row[j];
                }
            }
        }
        let norm = q.iter().map(|row| row[i] * row[i]).sum::<f64>().sqrt();
        if !(norm > 0.0 && norm.is_finite()) {
            return None;
        }
        r[i][i] = norm;
        for row in q.iter_mut() {
            row[i] /= norm;
        }
    }
    // B = Rᵀ Qᵀ: Rᵀ y = v (forward substitution), ξ = Q y.
    let y0 = v[0] / r[0][0];
    let y1 = (v[1] - r[0][1] * y0) / r[1][1];
    let y2 = (v[2] - r[0][2] * y0 - r[1][2] * y1) / r[2][2];
    let xi: Vec<f64> = q
        .iter()
        .map(|row| row[0] * y0 + row[1] * y1 + row[2] * y2)
        .collect();
    if !xi.iter().all(|x| x.is_finite()) {
        return None;
    }
    let mut u = [0.0; 3];
    let mut back = [0.0; 3];
    for (m, &x) in xi.iter().enumerate() {
        for c in 0..3 {
            u[c] += x * on_input[m][c];
            back[c] += x * on_triangulation[m][c];
        }
    }
    let residual = (0..3).map(|c| (back[c] - v[c]).abs()).fold(0.0, f64::max);
    Some((u, residual))
}

/// The result of App. A for a point of a region, before anything is done to coordinates
/// within rounding of zero.
#[derive(Debug, Clone, PartialEq)]
pub struct RecoveredPoint {
    /// The input face that contains the region.
    pub input_face: u32,
    /// The region.
    pub region: Region,
    /// The barycentric coordinates on the input face, as Eq. 16 gives them (they sum to one up
    /// to rounding; one may be slightly negative near an input edge).
    pub bary: [f64; 3],
    /// The number `ρ` of corners of the region used (all of them, 3 to 6, SPEC I16).
    pub n_corners: usize,
    /// The corners of the region on the input face (barycentric coordinates there).
    pub corners: Vec<[f64; 3]>,
    /// The residual of Eq. 15 (see [`recover_barycentric`]).
    pub residual: f64,
}

/// How the place of an inserted vertex was made a point strictly inside its input face
/// ([`RecoveredPoint::inside_input_face`], SPEC I16).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecoverReport {
    /// The input face.
    pub input_face: u32,
    /// The number of corners of the region used.
    pub n_corners: usize,
    /// The smallest barycentric coordinate as Eq. 16 gave it, before any move.
    pub min_coordinate: f64,
    /// The residual of Eq. 15.
    pub residual: f64,
    /// The fraction `λ` of the way to the centroid of the region's corners by which the point
    /// was moved to put it strictly inside the input face; `0` when it was inside.
    pub moved: f64,
}

impl RecoveredPoint {
    /// The smallest coordinate.
    pub fn min_coordinate(&self) -> f64 {
        self.bary[0].min(self.bary[1]).min(self.bary[2])
    }

    /// The place of a vertex inserted at this point (SPEC §2.2, I16): a point strictly inside
    /// the input face. Coordinates within [`RECOVER_ZERO`] of `(0, 1)` are rounding: the point
    /// is moved towards the centroid of the region's corners by the smallest fraction that puts
    /// every coordinate at or above `RECOVER_ZERO` times the centroid's (reported). `f` names
    /// the face of the triangulation in errors.
    ///
    /// # Errors
    /// [`IntegerError::RecoverFailed`] for a point farther outside, or when the move does not
    /// give a point strictly inside.
    pub fn inside_input_face(
        &self,
        f: usize,
    ) -> Result<(SurfacePoint, RecoverReport), IntegerError> {
        let w = self.bary;
        let min = self.min_coordinate();
        let max = w[0].max(w[1]).max(w[2]);
        let inside = |b: &[f64; 3]| b.iter().all(|&x| x > 0.0 && x < 1.0);
        let fail = |detail: String| IntegerError::RecoverFailed {
            face: f as u32,
            input_face: self.input_face,
            detail,
        };
        let mut moved = 0.0;
        let bary = if inside(&w) {
            w
        } else if min >= -RECOVER_ZERO && max <= 1.0 + RECOVER_ZERO {
            let rho = self.corners.len() as f64;
            let centroid = [0, 1, 2].map(|c| self.corners.iter().map(|x| x[c]).sum::<f64>() / rho);
            let mut lambda = 0.0f64;
            for c in 0..3 {
                let target = RECOVER_ZERO * centroid[c];
                if w[c] < target && centroid[c] > w[c] {
                    lambda = lambda.max((target - w[c]) / (centroid[c] - w[c]));
                }
            }
            let lambda = lambda.min(1.0);
            let b = [0, 1, 2].map(|c| (1.0 - lambda) * w[c] + lambda * centroid[c]);
            if !inside(&b) {
                return Err(fail(format!(
                    "the coordinates {w:?} on input face {} are within rounding of its boundary, \
                     and moving them towards the centroid {centroid:?} of the region's corners \
                     gives {b:?}",
                    self.input_face
                )));
            }
            moved = lambda;
            b
        } else {
            return Err(fail(format!(
                "the coordinates {w:?} on input face {} are outside it by more than {RECOVER_ZERO}",
                self.input_face
            )));
        };
        Ok((
            SurfacePoint::Face {
                face: self.input_face,
                bary,
            },
            RecoverReport {
                input_face: self.input_face,
                n_corners: self.n_corners,
                min_coordinate: min,
                residual: self.residual,
                moved,
            },
        ))
    }

    /// The point on the input mesh for a point query (SPEC §4.9): the point of the input face,
    /// or, when coordinates are within [`RECOVER_ZERO`] of zero, the point of the input edge or
    /// the input vertex they are within rounding of.
    ///
    /// # Errors
    /// [`IntegerError::RecoverFailed`] for a point farther outside.
    pub fn on_input(
        &self,
        tri: &IntegerTriangulation,
        f: usize,
    ) -> Result<SurfacePoint, IntegerError> {
        classify_in_input_face(tri, self.input_face as usize, self.bary).ok_or_else(|| {
            IntegerError::RecoverFailed {
                face: f as u32,
                input_face: self.input_face,
                detail: format!(
                    "the coordinates {:?} on input face {} are outside it by more than \
                     {RECOVER_ZERO}",
                    self.bary, self.input_face
                ),
            }
        })
    }
}

/// The point of input face `face` with coordinates `w`: a point of the face if all are above
/// [`RECOVER_ZERO`]; coordinates in `[−RECOVER_ZERO, RECOVER_ZERO]` are taken as zero (the rule
/// of WP5's `T⁰ → T¹` queries, `LOCATE_SNAP`), making it a point of an edge (one zero) or a
/// vertex (two); `None` farther outside.
fn classify_in_input_face(
    tri: &IntegerTriangulation,
    face: usize,
    w: [f64; 3],
) -> Option<SurfacePoint> {
    if w.iter().any(|x| !x.is_finite() || *x < -RECOVER_ZERO) {
        return None;
    }
    let mesh = tri.input.mesh();
    let zero = [0, 1, 2].map(|c| w[c] <= RECOVER_ZERO);
    match zero.iter().filter(|&&z| z).count() {
        0 => {
            let s = w[0] + w[1] + w[2];
            Some(SurfacePoint::Face {
                face: face as u32,
                bary: [w[0] / s, w[1] / s, w[2] / s],
            })
        }
        1 => {
            // On the side opposite corner k: halfedge k + 1, from corner k + 1 to k + 2.
            let k = zero.iter().position(|&z| z)?;
            let (i, j) = ((k + 1) % 3, (k + 2) % 3);
            let h = 3 * face + i;
            let t_h = w[j] / (w[i] + w[j]);
            let edge = mesh.edge(h);
            let forward = mesh.from_vertex(h) == mesh.edge_vertices(edge)[0];
            Some(SurfacePoint::Edge {
                edge: edge as u32,
                t: if forward { t_h } else { 1.0 - t_h },
            })
        }
        _ => {
            let k = (0..3).max_by(|&a, &b| w[a].total_cmp(&w[b]))?;
            Some(SurfacePoint::Vertex(mesh.faces()[face][k]))
        }
    }
}

/// The barycentric coordinates on input face `face` of a point of the input mesh, if it is a
/// point of the face's closure.
pub(crate) fn bary_in_input_face(
    tri: &IntegerTriangulation,
    face: usize,
    q: SurfacePoint,
) -> Option<[f64; 3]> {
    let mesh = tri.input.mesh();
    let corners = mesh.faces()[face];
    let mut b = [0.0; 3];
    match q {
        SurfacePoint::Vertex(v) => b[corners.iter().position(|&c| c == v)?] = 1.0,
        SurfacePoint::Edge { edge, t } => {
            let s = (0..3).find(|&s| mesh.edge(3 * face + s) == edge as usize)?;
            let forward = mesh.from_vertex(3 * face + s) == mesh.edge_vertices(edge as usize)[0];
            let th = if forward { t } else { 1.0 - t };
            b[s] = 1.0 - th;
            b[(s + 1) % 3] = th;
        }
        SurfacePoint::Face { face: g, bary } => {
            if g as usize != face {
                return None;
            }
            b = bary;
        }
    }
    Some(b)
}

/// App. A for the point `v` of `region` of the face `arcs`: the input face of the region, its
/// corners there, Eq. 15–16.
pub(crate) fn recover_in_input_face(
    tri: &IntegerTriangulation,
    arcs: &FaceArcs,
    region: Region,
    v: [f64; 3],
) -> Result<RecoveredPoint, IntegerError> {
    let f = arcs.face;
    let input_face = tri.region_input_face(arcs, region)?;
    let points = arcs.region_points(region);
    let fail = |detail: String| IntegerError::RecoverFailed {
        face: f as u32,
        input_face: input_face as u32,
        detail,
    };
    let mut on_t1 = Vec::with_capacity(points.len());
    let mut on_t0 = Vec::with_capacity(points.len());
    for &p in &points {
        let q = match p {
            RegionPoint::Corner(k) => tri.locations[tri.complex.faces()[f][k] as usize],
            RegionPoint::Crossing(k, i) => arcs.sides[k].points[i],
        };
        let b = bary_in_input_face(tri, input_face, q).ok_or_else(|| {
            fail(format!(
                "the corner {p:?} of region {region} is at {q}, not on input face {input_face}"
            ))
        })?;
        on_t1.push(arcs.bary_of(p));
        on_t0.push(b);
    }
    let (bary, residual) = recover_barycentric(&on_t1, &on_t0, v).ok_or_else(|| {
        fail(format!(
            "the {} corners of region {region} do not span the face",
            points.len()
        ))
    })?;
    Ok(RecoveredPoint {
        input_face: input_face as u32,
        region,
        bary,
        n_corners: points.len(),
        corners: on_t0,
        residual,
    })
}

impl IntegerTriangulation {
    /// The input face that contains `region` of the face `arcs` (see the
    /// [module documentation](crate::recover)).
    fn region_input_face(&self, arcs: &FaceArcs, region: Region) -> Result<usize, IntegerError> {
        let f = arcs.face;
        let mesh = self.input.mesh();
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Face(f as u32), detail)
        };
        let Some((k, slot)) = arcs.region_gap(region) else {
            // All three sides run along input edges: the face is the input face on the left of
            // the input halfedge along a side that leaves an input vertex.
            for k in 0..3 {
                let h = 3 * f + k;
                let a = self.complex.from_vertex(h);
                if !self.is_input_vertex(a) {
                    continue;
                }
                let r = self.roundabouts[h];
                if r >= self.local.degree(a) {
                    return Err(invalid(format!(
                        "the roundabout {r} of halfedge {h} is out of range"
                    )));
                }
                let edge = self.local.edge(a, r);
                return mesh
                    .outgoing_halfedges(a)
                    .find(|&g| mesh.edge(g) == edge && self.local.of_halfedge(g) == r)
                    .map(|g| g / 3)
                    .ok_or_else(|| {
                        invalid(format!(
                            "side {k} runs along input edge {edge}, which has no halfedge \
                             leaving vertex {a}"
                        ))
                    });
            }
            return Err(invalid(
                "all three sides run along input edges, but no corner is an input vertex".into(),
            ));
        };
        let side = &arcs.sides[k];
        let trace = side
            .trace
            .as_ref()
            .ok_or_else(|| invalid(format!("side {k} has no path over the input mesh")))?;
        let x = &trace.crossings;
        let n = x.len();
        if slot > n {
            return Err(invalid(format!(
                "slot {slot} of side {k} with {n} crossings"
            )));
        }
        let q = if side.reversed { n - slot } else { slot };
        if q < n {
            Ok(x[q].halfedge as usize / 3)
        } else if n > 0 {
            mesh.twin(x[n - 1].halfedge as usize)
                .map(|g| g / 3)
                .ok_or_else(|| {
                    invalid(format!(
                        "side {k} leaves the input mesh at its last crossing"
                    ))
                })
        } else {
            trace
                .first_face
                .map(|g| g as usize)
                .ok_or_else(|| invalid(format!("side {k} has no crossing and no first input face")))
        }
    }

    /// The place on the input mesh of a point of the triangulation (SPEC §4.9; the signature of
    /// [`ddg_intrinsic::IntrinsicTriangulation::point_on_input`]):
    ///
    /// * a vertex: its place;
    /// * a point of an edge along an input edge (`n = −1`): interpolated between the places of
    ///   its ends on that input edge;
    /// * a point of an edge that crosses input edges: between the two neighbouring crossings
    ///   (from `ExtractIntrinsicEdge`), interpolated in their input face;
    /// * a point inside a face: its region and App. A ([`locate_in_face`](Self::locate_in_face)).
    ///
    /// A point of a face with a zero coordinate is the point of the side it lies on; one with a
    /// coordinate `1` is that corner. A result within [`RECOVER_ZERO`] of a side or corner of its
    /// input face is returned as that input edge's (or vertex's) point.
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`ddg_intrinsic::IntrinsicError::InvalidPoint`] for a point that is not one of the
    /// triangulation; [`IntegerError::RecoverFailed`]; the errors of
    /// [`extract_intrinsic_edge`](Self::extract_intrinsic_edge) and
    /// [`locate_in_face`](Self::locate_in_face).
    pub fn point_on_input(&self, point: SurfacePoint) -> Result<SurfacePoint, IntegerError> {
        let invalid_point = |reason: &'static str| {
            IntegerError::Intrinsic(ddg_intrinsic::IntrinsicError::InvalidPoint { point, reason })
        };
        let c = &self.complex;
        match point {
            SurfacePoint::Vertex(v) => {
                let nv = c.n_vertices();
                if v as usize >= nv {
                    return Err(IntegerError::IndexOutOfRange {
                        what: "vertex",
                        index: v as usize,
                        bound: nv,
                    });
                }
                Ok(self.locations[v as usize])
            }
            SurfacePoint::Edge { edge, t } => {
                let e = edge as usize;
                if e >= c.n_edges() {
                    return Err(IntegerError::IndexOutOfRange {
                        what: "edge",
                        index: e,
                        bound: c.n_edges(),
                    });
                }
                if !(t.is_finite() && (0.0..=1.0).contains(&t)) {
                    return Err(invalid_point("the edge parameter is not in [0, 1]"));
                }
                let h = c.edge_halfedge(e);
                if t == 0.0 {
                    return Ok(self.locations[c.from_vertex(h)]);
                }
                if t == 1.0 {
                    return Ok(self.locations[c.to_vertex(h)]);
                }
                if self.normal[e] < 0 {
                    self.point_along_input_edge(e, t)
                } else {
                    self.point_on_crossed_edge(e, t)
                }
            }
            SurfacePoint::Face { face, bary } => {
                let f = face as usize;
                if f >= c.n_faces() {
                    return Err(IntegerError::IndexOutOfRange {
                        what: "face",
                        index: f,
                        bound: c.n_faces(),
                    });
                }
                if !bary.iter().all(|b| b.is_finite() && *b >= 0.0)
                    || ((bary[0] + bary[1] + bary[2]) - 1.0).abs() > 1e-9
                {
                    return Err(invalid_point(
                        "the barycentric coordinates are not non-negative numbers summing to one",
                    ));
                }
                if let Some(k) = bary.iter().position(|&b| b == 1.0) {
                    return Ok(self.locations[c.faces()[f][k] as usize]);
                }
                if let Some(k) = bary.iter().position(|&b| b == 0.0) {
                    // On the side opposite corner k: halfedge k + 1.
                    let h = 3 * f + (k + 1) % 3;
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    let th = bary[j] / (bary[i] + bary[j]);
                    let e = c.edge(h);
                    let t = if c.edge_halfedge(e) == h {
                        th
                    } else {
                        1.0 - th
                    };
                    return self.point_on_input(SurfacePoint::Edge { edge: e as u32, t });
                }
                let location = self.locate_in_face(f, bary)?;
                location.recovered.on_input(self, f)
            }
        }
    }

    /// The input edge along which edge `e` (`n = −1`) runs: named by the place of an inserted
    /// end, else by the roundabout at an input-vertex end.
    fn input_edge_of_along(&self, e: usize) -> Result<usize, IntegerError> {
        let c = &self.complex;
        let h = c.edge_halfedge(e);
        let (a, b) = (c.from_vertex(h), c.to_vertex(h));
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::AlongInputEdge, Place::Edge(e as u32), detail)
        };
        for v in [a, b] {
            if let SurfacePoint::Edge { edge, .. } = self.locations[v] {
                return Ok(edge as usize);
            }
        }
        for (v, g) in [(a, Some(h)), (b, c.twin(h))] {
            if let (true, Some(g)) = (self.is_input_vertex(v), g) {
                let r = self.roundabouts[g];
                if r < self.local.degree(v) {
                    return Ok(self.local.edge(v, r));
                }
                return Err(invalid(format!(
                    "the roundabout {r} of halfedge {g} is out of range"
                )));
            }
        }
        Err(invalid("no end of the edge names an input edge".into()))
    }

    /// The point at parameter `t` of edge `e` along an input edge (`n = −1`).
    fn point_along_input_edge(&self, e: usize, t: f64) -> Result<SurfacePoint, IntegerError> {
        let c = &self.complex;
        let h = c.edge_halfedge(e);
        let input_edge = self.input_edge_of_along(e)?;
        let ta = self.parameter_on_input_edge(c.from_vertex(h), input_edge)?;
        let tb = self.parameter_on_input_edge(c.to_vertex(h), input_edge)?;
        let s = (1.0 - t) * ta + t * tb;
        let [lo, hi] = self.input.mesh().edge_vertices(input_edge);
        Ok(if s <= 0.0 {
            SurfacePoint::Vertex(lo as u32)
        } else if s >= 1.0 {
            SurfacePoint::Vertex(hi as u32)
        } else {
            SurfacePoint::Edge {
                edge: input_edge as u32,
                t: s,
            }
        })
    }

    /// The point at parameter `t` of edge `e` (`n ≥ 0`): between the neighbouring crossings, in
    /// their input face.
    fn point_on_crossed_edge(&self, e: usize, t: f64) -> Result<SurfacePoint, IntegerError> {
        let trace = self.extract_intrinsic_edge(e, &EdgeTraceOptions::default())?;
        let mesh = self.input.mesh();
        let invalid = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Edge(e as u32), detail)
        };
        let at = t * self.lengths[e];
        let x = &trace.crossings;
        let before = x.iter().filter(|c| c.distance < at).count();
        let face = match (x.get(before), before.checked_sub(1).map(|k| &x[k])) {
            (Some(next), _) => next.halfedge as usize / 3,
            (None, Some(last)) => {
                mesh.twin(last.halfedge as usize)
                    .ok_or_else(|| invalid("its last crossing is a boundary input edge".into()))?
                    / 3
            }
            (None, None) => trace
                .first_face
                .ok_or_else(|| invalid("it has no crossing and no first input face".into()))?
                as usize,
        };
        let h = self.complex.edge_halfedge(e);
        let crossing = |k: usize| SurfacePoint::Edge {
            edge: x[k].edge,
            t: x[k].t,
        };
        let (p0, d0) = match before.checked_sub(1) {
            Some(k) => (crossing(k), x[k].distance),
            None => (self.locations[self.complex.from_vertex(h)], 0.0),
        };
        let (p1, d1) = match x.get(before) {
            Some(c) => (crossing(before), c.distance),
            None => (self.locations[self.complex.to_vertex(h)], self.lengths[e]),
        };
        let (Some(b0), Some(b1)) = (
            bary_in_input_face(self, face, p0),
            bary_in_input_face(self, face, p1),
        ) else {
            return Err(invalid(format!(
                "the points {p0} and {p1} around the point are not on input face {face}"
            )));
        };
        let s = if d1 > d0 {
            ((at - d0) / (d1 - d0)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let w = [0, 1, 2].map(|k| (1.0 - s) * b0[k] + s * b1[k]);
        classify_in_input_face(self, face, w).ok_or_else(|| {
            invalid(format!(
                "the interpolated point {w:?} is outside input face {face}"
            ))
        })
    }
}
