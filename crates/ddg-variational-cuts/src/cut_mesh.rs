//! The cut mesh (SPEC §3.3): the input mesh with every crossed face split along its cut segment,
//! and the two sides of the cut separated. [`CutTriangles`] is the raw form the evaluation works
//! with, [`CutMesh`] the same mesh as a validated mesh with geometry.

use ddg_core::{Geometry, GeometryError, TriMesh};

use crate::cut::{Cut, lone_corner};
use crate::error::CutError;
use crate::problem::{UnitAreaScale, check_length};

/// The inner product of two tangent vectors of a triangle `abc` given as differences of
/// barycentric coordinates `d = (d_a, d_b, d_c)`, `e = (e_a, e_b, e_c)` (each summing to zero),
/// from the squared edge lengths `len2 = [ℓ_ab², ℓ_bc², ℓ_ca²]`:
///
/// ```text
/// ⟨d, e⟩ = −½ [ (d_a e_b + d_b e_a) ℓ_ab² + (d_b e_c + d_c e_b) ℓ_bc² + (d_c e_a + d_a e_c) ℓ_ca² ]
/// ```
///
/// With `e = d` this is the squared distance of SPEC §3.3,
/// `|x − y|² = −(d_a d_b ℓ_ab² + d_b d_c ℓ_bc² + d_c d_a ℓ_ca²)`.
///
/// [formula: cut-barycentric-length]
#[inline]
pub fn barycentric_dot(d: [f64; 3], e: [f64; 3], len2: [f64; 3]) -> f64 {
    -0.5 * ((d[0] * e[1] + d[1] * e[0]) * len2[0]
        + (d[1] * e[2] + d[2] * e[1]) * len2[1]
        + (d[2] * e[0] + d[0] * e[2]) * len2[2])
}

/// The area of the triangle whose corners have barycentric coordinates `x[0]`, `x[1]`, `x[2]`
/// in a face of area `face_area`: `det(x) · face_area`. It is positive when the corners are in
/// the orientation of the face, and the areas of a subdivision of the face sum to `face_area`.
///
/// [formula: cut-subtriangle-area]
#[inline]
pub fn subtriangle_area(x: [[f64; 3]; 3], face_area: f64) -> f64 {
    let det = x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1])
        - x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0])
        + x[0][2] * (x[1][0] * x[2][1] - x[1][1] * x[2][0]);
    det * face_area
}

/// The aspect ratio of SPEC §3.7 from the squared edge lengths and the area of a triangle:
/// `(longest edge)² · √3 / (4 · area)`, which is 1 for an equilateral triangle.
#[inline]
pub fn aspect_ratio(len2: [f64; 3], area: f64) -> f64 {
    len2[0].max(len2[1]).max(len2[2]) * 3f64.sqrt() / (4.0 * area)
}

/// Per-face data of the input mesh on the unit-area surface.
#[derive(Debug, Clone)]
pub(crate) struct BaseGeometry {
    /// Squared length of halfedge `3 f + k`.
    pub len2: Vec<[f64; 3]>,
    /// Face areas; they sum to 1.
    pub area: Vec<f64>,
    /// Corner cotangents, laid out like `Geometry::corner_cotans`.
    pub cot: Vec<[f64; 3]>,
}

impl BaseGeometry {
    pub(crate) fn new(mesh: &TriMesh, geometry: &Geometry, scale: &UnitAreaScale) -> BaseGeometry {
        let lengths = geometry.edge_lengths();
        let len2 = (0..mesh.n_faces())
            .map(|f| {
                [0, 1, 2].map(|k| {
                    let l = lengths[mesh.edge(3 * f + k)] * scale.length_scale;
                    l * l
                })
            })
            .collect();
        let area = geometry
            .face_areas()
            .iter()
            .map(|a| a / scale.total_area)
            .collect();
        BaseGeometry {
            len2,
            area,
            cot: geometry.corner_cotans().to_vec(),
        }
    }
}

/// What an edge of a cut-mesh triangle is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum EdgeKind {
    /// Shared with another triangle of the same region.
    Interior,
    /// A cut segment: a boundary edge of the region, on the cut.
    Cut,
    /// A boundary edge of the input mesh, or the part of one on this side of a cut point.
    Boundary,
}

/// Where a vertex of the cut mesh comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutVertex {
    /// A vertex of the input mesh.
    Original(u32),
    /// One side of a node of the cut: a cut point (an index into [`Cut::points`]) or, from
    /// `Cut::points.len()` on, a triple point (SPEC §8.2; never for two regions).
    Point(u32),
}

/// The cut mesh (SPEC §3.3) as a list of triangles with their intrinsic geometry on the
/// unit-area surface: the input mesh with every crossed face split along its cut segment, and
/// the two sides of the cut separated.
///
/// # Vertices
///
/// The cut mesh is the cut *surface*: a cut point has one vertex on each side of the cut, so
/// every vertex of the cut mesh lies in exactly one region. With `n_V` vertices in the input
/// mesh,
///
/// * vertex `i < n_V` is vertex `i` of the input mesh;
/// * cut point `k` on the edge `[a, b]` (`a < b`) has the two copies `n_V + 2k` (in the region
///   of `a`) and `n_V + 2k + 1` (in the region of `b`);
/// * triple point `m` in the face with the corners `(c_0, c_1, c_2)` has, after all cut points,
///   the three copies `n_V + 2 n_P + 3m + j` in the region of `c_j` (SPEC §8.2).
///
/// Use [`CutTriangles::vertex_source`], [`CutTriangles::point_vertices`] and
/// [`CutTriangles::triple_vertices`] instead of relying on this layout.
///
/// # Triangles
///
/// Triangles are listed in order of their parent face: an uncut face gives one triangle (the
/// face itself, same corners), a face with three regions six (for every corner `c` in turn,
/// with `p` the cut point on the edge that leaves `c`, `p′` the one on the edge that arrives
/// at `c` and `m` the triple point: `(c, p, m)`, then `(c, m, p′)`), a face with two regions
/// three: with `a` the vertex that is alone on its side
/// and `p ∈ ab`, `q ∈ ca` the cut points, `(a, p, q)`, then `(p, b, c)`, `(p, c, q)` if
/// `|pc| ≤ |bq|`, else `(p, b, q)`, `(b, c, q)`. Corner slot `k` of triangle `t` is the corner
/// at vertex `faces[t][k]`; edge slot `k` is the edge from corner `k` to corner `k + 1`.
/// Orientation is inherited from the parent face.
///
/// # Two forms
///
/// This is the raw form the evaluation works with: one record per triangle with its
/// cotangents, area and squared edge lengths, computed in the plane of the parent face from
/// barycentric coordinates ([`barycentric_dot`], [`subtriangle_area`]). It is cheap to build.
/// [`CutMesh`] is the same mesh as a validated [`TriMesh`] with a [`Geometry`], which is what a
/// flattening consumes; it is built on demand ([`CutTriangles::to_mesh`]).
#[derive(Debug, Clone, PartialEq)]
pub struct CutTriangles {
    /// Number of vertices of the input mesh.
    pub n_original: usize,
    /// Number of cut points.
    pub n_points: usize,
    /// Number of triple points (SPEC §8.2); zero for two regions.
    pub n_triples: usize,
    /// The triangles, as vertex ids of the cut mesh.
    pub faces: Vec<[u32; 3]>,
    /// The face of the input mesh each triangle is part of.
    pub parent: Vec<u32>,
    /// The region of each triangle.
    pub region: Vec<u32>,
    /// The cotangent of each corner angle (scale invariant).
    pub cot: Vec<[f64; 3]>,
    /// The area of each triangle, on the unit-area surface.
    pub area: Vec<f64>,
    /// The squared length of each edge (edge `k` runs from corner `k` to corner `k + 1`), on the
    /// unit-area surface.
    pub len2: Vec<[f64; 3]>,
    /// What each edge is (edge slots as for `len2`).
    pub edge_kind: Vec<[EdgeKind; 3]>,
    /// The region of every vertex of the cut mesh.
    pub vertex_region: Vec<u32>,
    /// The length of every cut segment (aligned with [`Cut::segments`]), on the unit-area
    /// surface.
    pub segment_length: Vec<f64>,
    /// For every cut segment, the triangle `(a, p, q)` that it cuts off its face: `a` is the
    /// vertex of the face that is alone on its side, and the segment is edge slot 1. For a
    /// segment from the cut point `p` to a triple point `m`: the triangle `(c, p, m)` with `c`
    /// the corner that the edge of `p` leaves; the segment is edge slot 1 again.
    pub segment_triangle: Vec<u32>,
    /// For every cut point, whether an incident triangle of the cut mesh has an aspect ratio
    /// above the threshold it was built with (SPEC §3.7).
    pub sliver_point: Vec<bool>,
    /// The largest aspect ratio among the triangles of cut faces (1 if no face is cut).
    pub max_aspect_ratio: f64,
    /// Number of cut faces in which the squared lengths of the two diagonals of the quad agree
    /// to a relative `10⁻⁹`, so that rounding decides which one splits it (SPEC §3.3). The result of the evaluation
    /// depends on that choice only at the level of the discretization error.
    pub n_diagonal_near_ties: usize,
}

/// The cut mesh as a validated mesh with geometry and maps back to the input (SPEC §7): what a
/// flattening consumes. Vertices and faces are numbered as in the [`CutTriangles`] it was built
/// from.
#[derive(Debug, Clone, PartialEq)]
pub struct CutMesh {
    /// Connectivity. Its boundary is the cut (twice, once per side) and the boundary of the
    /// input mesh.
    pub mesh: TriMesh,
    /// Edge lengths, angles and areas, in the length unit of the caller's geometry (not the
    /// unit-area surface).
    pub geometry: Geometry,
    /// For every face, the face of the input mesh it is part of.
    pub face_parent: Vec<u32>,
    /// For every face, its region.
    pub face_region: Vec<u32>,
    /// For every vertex, the vertex of the input mesh or the cut point it is.
    pub vertex_source: Vec<CutVertex>,
    /// For every vertex, its region.
    pub vertex_region: Vec<u32>,
}

/// Geometry of the triangle with barycentric corners `x` in a face with squared edge lengths
/// `len2 = [ℓ_ab², ℓ_bc², ℓ_ca²]` and area `face_area`: `(area, squared edge lengths, corner
/// cotangents)`, or `None` if it has no positive area.
#[inline]
fn subtriangle(
    x: [[f64; 3]; 3],
    len2: [f64; 3],
    face_area: f64,
) -> Option<(f64, [f64; 3], [f64; 3])> {
    let area = subtriangle_area(x, face_area);
    // d[k] is edge k, from corner k to corner k + 1.
    let d = [0, 1, 2].map(|k| {
        let (p, q) = (x[k], x[(k + 1) % 3]);
        [q[0] - p[0], q[1] - p[1], q[2] - p[2]]
    });
    let sub_len2 = d.map(|e| barycentric_dot(e, e, len2));
    // cot at corner k = <edge k, −edge k−1> / (2 area).
    let cot = [0, 1, 2].map(|k| -barycentric_dot(d[k], d[(k + 2) % 3], len2) / (2.0 * area));
    let valid = area > 0.0
        && area.is_finite()
        && sub_len2.iter().all(|l| *l > 0.0 && l.is_finite())
        && cot.iter().all(|c| c.is_finite());
    valid.then_some((area, sub_len2, cot))
}

impl CutTriangles {
    /// Builds the cut mesh of `cut` on the unit-area geometry `base` (SPEC §3.3): a crossed face
    /// `abc` with `a` alone on its side and cut points `p ∈ ab`, `q ∈ ca` becomes the triangle
    /// `(a, p, q)` and the quad `(p, b, c, q)`, which is split by the shorter of its diagonals
    /// `p–c`, `b–q` (tie: `p–c`).
    pub(crate) fn build(
        mesh: &TriMesh,
        base: &BaseGeometry,
        cut: &Cut,
        sliver_aspect_ratio: f64,
    ) -> Result<CutTriangles, CutError> {
        let nv = mesh.n_vertices();
        let n_points = cut.points.len();
        let n_triangles = mesh.n_faces() + 2 * cut.segments.len();
        let mut next_triple = 0usize;
        let mut out = CutTriangles {
            n_original: nv,
            n_points,
            n_triples: cut.triple_points.len(),
            faces: Vec::with_capacity(n_triangles),
            parent: Vec::with_capacity(n_triangles),
            region: Vec::with_capacity(n_triangles),
            cot: Vec::with_capacity(n_triangles),
            area: Vec::with_capacity(n_triangles),
            len2: Vec::with_capacity(n_triangles),
            edge_kind: Vec::with_capacity(n_triangles),
            vertex_region: Vec::with_capacity(nv + 2 * n_points),
            segment_length: Vec::with_capacity(cut.segments.len()),
            segment_triangle: Vec::with_capacity(cut.segments.len()),
            sliver_point: vec![false; n_points],
            max_aspect_ratio: 1.0,
            n_diagonal_near_ties: 0,
        };
        out.vertex_region.extend_from_slice(&cut.vertex_region);
        for point in &cut.points {
            let [a, b] = mesh.edge_vertices(point.edge as usize);
            out.vertex_region.push(cut.vertex_region[a]);
            out.vertex_region.push(cut.vertex_region[b]);
        }
        for triple in &cut.triple_points {
            for v in mesh.faces()[triple.face as usize] {
                out.vertex_region.push(cut.vertex_region[v as usize]);
            }
        }

        let boundary_kind = |h: usize| {
            if mesh.twin(h).is_some() {
                EdgeKind::Interior
            } else {
                EdgeKind::Boundary
            }
        };
        for (f, face) in mesh.faces().iter().enumerate() {
            let regions = face.map(|v| cut.vertex_region[v as usize]);
            let k = match lone_corner(regions) {
                Ok(None) => {
                    out.faces.push(*face);
                    out.parent.push(f as u32);
                    out.region.push(regions[0]);
                    out.cot.push(base.cot[f]);
                    out.area.push(base.area[f]);
                    out.len2.push(base.len2[f]);
                    out.edge_kind
                        .push([0, 1, 2].map(|k| boundary_kind(3 * f + k)));
                    continue;
                }
                Ok(Some(k)) => k,
                Err(()) => {
                    // Three regions (SPEC §8.2): every corner keeps the quad between its two cut
                    // points and the triple point, split by the diagonal to the triple point.
                    let triple = next_triple;
                    next_triple += 1;
                    let len2 = base.len2[f];
                    let xm = cut.triple_points[triple].barycentric;
                    let mut on_edge = [(0u32, [0.0; 3]); 3];
                    for (k, slot) in on_edge.iter_mut().enumerate() {
                        let e = mesh.edge(3 * f + k);
                        let point = cut.edge_point[e];
                        let t = cut.points[point as usize].t;
                        let s = if mesh.edges()[e][0] == face[k] {
                            t
                        } else {
                            1.0 - t
                        };
                        let mut x = [0.0; 3];
                        x[k] = 1.0 - s;
                        x[(k + 1) % 3] = s;
                        *slot = (point, x);
                    }
                    for c in 0..3 {
                        let previous = (c + 2) % 3;
                        let (p, xp) = on_edge[c];
                        let (q, xq) = on_edge[previous];
                        let copy = |point: u32| -> u32 {
                            let lo = mesh.edges()[cut.points[point as usize].edge as usize][0];
                            (nv + 2 * point as usize + usize::from(face[c] != lo)) as u32
                        };
                        let m = (nv + 2 * n_points + 3 * triple + c) as u32;
                        let mut xc = [0.0; 3];
                        xc[c] = 1.0;
                        let parts = [
                            (
                                [face[c], copy(p), m],
                                [xc, xp, xm],
                                [boundary_kind(3 * f + c), EdgeKind::Cut, EdgeKind::Interior],
                                p,
                            ),
                            (
                                [face[c], m, copy(q)],
                                [xc, xm, xq],
                                [
                                    EdgeKind::Interior,
                                    EdgeKind::Cut,
                                    boundary_kind(3 * f + previous),
                                ],
                                q,
                            ),
                        ];
                        for (part, (ids, x, kinds, point)) in parts.into_iter().enumerate() {
                            let Some((area, sub_len2, cot)) = subtriangle(x, len2, base.area[f])
                            else {
                                return Err(CutError::DegenerateCutTriangle {
                                    face: f as u32,
                                    vertices: *face,
                                    detail: format!(
                                        "three regions meet in this face, and the triangle at \
                                         vertex {} between a cut point and the triple point has \
                                         no area",
                                        face[c]
                                    ),
                                });
                            };
                            let aspect = aspect_ratio(sub_len2, area);
                            out.max_aspect_ratio = out.max_aspect_ratio.max(aspect);
                            if aspect > sliver_aspect_ratio {
                                out.sliver_point[point as usize] = true;
                            }
                            if part == 0 {
                                // The segment from the cut point on the edge leaving c to m.
                                out.segment_triangle.push(out.faces.len() as u32);
                                out.segment_length.push(sub_len2[1].sqrt());
                            }
                            out.faces.push(ids);
                            out.parent.push(f as u32);
                            out.region.push(regions[c]);
                            out.cot.push(cot);
                            out.area.push(area);
                            out.len2.push(sub_len2);
                            out.edge_kind.push(kinds);
                        }
                    }
                    continue;
                }
            };
            // Rotate the face so that the lone vertex comes first: (a, b, c), with halfedges
            // a -> b, b -> c, c -> a.
            let (ka, kb, kc) = (k, (k + 1) % 3, (k + 2) % 3);
            let (a, b, c) = (face[ka], face[kb], face[kc]);
            let (h_ab, h_bc, h_ca) = (3 * f + ka, 3 * f + kb, 3 * f + kc);
            let len2 = [base.len2[f][ka], base.len2[f][kb], base.len2[f][kc]];
            let face_area = base.area[f];
            let (kind_ab, kind_bc, kind_ca) = (
                boundary_kind(h_ab),
                boundary_kind(h_bc),
                boundary_kind(h_ca),
            );
            // The cut points and their parameters measured from a.
            let from_a = |h: usize| -> (u32, f64) {
                let e = mesh.edge(h);
                let point = cut.edge_point[e];
                let t = cut.points[point as usize].t;
                let s = if mesh.edges()[e][0] == a { t } else { 1.0 - t };
                (point, s)
            };
            let (p, sp) = from_a(h_ab);
            let (q, sq) = from_a(h_ca);
            let copy = |point: u32, side: u32| -> u32 {
                let lo = mesh.edges()[cut.points[point as usize].edge as usize][0];
                (nv + 2 * point as usize + usize::from(side != lo)) as u32
            };
            let (xa, xb, xc) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
            let xp = [1.0 - sp, sp, 0.0];
            let xq = [1.0 - sq, 0.0, sq];

            // The corner triangle (a, p, q) comes first.
            out.segment_triangle.push(out.faces.len() as u32);
            let mut push = |ids: [u32; 3],
                            x: [[f64; 3]; 3],
                            region: u32,
                            kinds: [EdgeKind; 3],
                            points: [Option<u32>; 2]|
             -> Result<[f64; 3], CutError> {
                let Some((area, sub_len2, cot)) = subtriangle(x, len2, face_area) else {
                    return Err(CutError::DegenerateCutTriangle {
                        face: f as u32,
                        vertices: *face,
                        detail: format!(
                            "the cut crosses the edges ({a}, {b}) and ({c}, {a}) at the fractions \
                             {sp} and {sq} of their length from vertex {a}"
                        ),
                    });
                };
                let aspect = aspect_ratio(sub_len2, area);
                out.max_aspect_ratio = out.max_aspect_ratio.max(aspect);
                if aspect > sliver_aspect_ratio {
                    for point in points.into_iter().flatten() {
                        out.sliver_point[point as usize] = true;
                    }
                }
                out.faces.push(ids);
                out.parent.push(f as u32);
                out.region.push(region);
                out.cot.push(cot);
                out.area.push(area);
                out.len2.push(sub_len2);
                out.edge_kind.push(kinds);
                Ok(sub_len2)
            };

            let (region_a, region_bc) = (regions[ka], regions[kb]);
            let (pa, qa) = (copy(p, a), copy(q, a));
            let (pb, qc) = (copy(p, b), copy(q, c));
            let corner = push(
                [a, pa, qa],
                [xa, xp, xq],
                region_a,
                [kind_ab, EdgeKind::Cut, kind_ca],
                [Some(p), Some(q)],
            )?;
            out.segment_length.push(corner[1].sqrt());
            let d_pc = [xc[0] - xp[0], xc[1] - xp[1], xc[2] - xp[2]];
            let d_bq = [xq[0] - xb[0], xq[1] - xb[1], xq[2] - xb[2]];
            let (pc2, bq2) = (
                barycentric_dot(d_pc, d_pc, len2),
                barycentric_dot(d_bq, d_bq, len2),
            );
            if (pc2 - bq2).abs() <= 1e-9 * pc2.max(bq2) {
                out.n_diagonal_near_ties += 1;
            }
            if pc2 <= bq2 {
                push(
                    [pb, b, c],
                    [xp, xb, xc],
                    region_bc,
                    [kind_ab, kind_bc, EdgeKind::Interior],
                    [Some(p), None],
                )?;
                push(
                    [pb, c, qc],
                    [xp, xc, xq],
                    region_bc,
                    [EdgeKind::Interior, kind_ca, EdgeKind::Cut],
                    [Some(p), Some(q)],
                )?;
            } else {
                push(
                    [pb, b, qc],
                    [xp, xb, xq],
                    region_bc,
                    [kind_ab, EdgeKind::Interior, EdgeKind::Cut],
                    [Some(p), Some(q)],
                )?;
                push(
                    [b, c, qc],
                    [xb, xc, xq],
                    region_bc,
                    [kind_bc, kind_ca, EdgeKind::Interior],
                    [None, Some(q)],
                )?;
            }
        }
        Ok(out)
    }

    /// Number of vertices of the cut mesh.
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.n_original + 2 * self.n_points + 3 * self.n_triples
    }

    /// Number of triangles of the cut mesh.
    #[inline]
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Where vertex `v` of the cut mesh comes from. A copy of triple point `m` is reported as
    /// `CutVertex::Point(n_points + m)`: the node numbering of
    /// [`CutSegment::points`](crate::CutSegment::points).
    #[inline]
    pub fn vertex_source(&self, v: usize) -> CutVertex {
        let first_triple = self.n_original + 2 * self.n_points;
        if v < self.n_original {
            CutVertex::Original(v as u32)
        } else if v < first_triple {
            CutVertex::Point(((v - self.n_original) / 2) as u32)
        } else {
            CutVertex::Point((self.n_points + (v - first_triple) / 3) as u32)
        }
    }

    /// The three vertices of the cut mesh at triple point `triple`: its copies in the regions of
    /// corner 0, 1 and 2 of its face.
    #[inline]
    pub fn triple_vertices(&self, triple: usize) -> [usize; 3] {
        let v = self.n_original + 2 * self.n_points + 3 * triple;
        [v, v + 1, v + 2]
    }

    /// The two vertices of the cut mesh at cut point `point`: the copy on the side of the
    /// lower end `a` of its edge `[a, b]`, then the copy on the side of `b`.
    #[inline]
    pub fn point_vertices(&self, point: usize) -> [usize; 2] {
        let v = self.n_original + 2 * point;
        [v, v + 1]
    }

    /// The two vertices of the cut mesh at cut point `point`, ordered by region: the copy in the
    /// region with the lower id first.
    #[inline]
    pub fn point_vertices_by_region(&self, point: usize) -> [usize; 2] {
        let [v0, v1] = self.point_vertices(point);
        if self.vertex_region[v0] <= self.vertex_region[v1] {
            [v0, v1]
        } else {
            [v1, v0]
        }
    }

    /// Extends a function on the vertices of the input mesh to the vertices of the cut mesh by
    /// linear interpolation along the crossed edges; both copies of a cut point get the same
    /// value.
    pub(crate) fn extend(&self, mesh: &TriMesh, cut: &Cut, values: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.n_vertices());
        out.extend_from_slice(values);
        for (node, value) in cut.interpolate_nodes(mesh, values).into_iter().enumerate() {
            let copies = if node < self.n_points { 2 } else { 3 };
            out.extend(std::iter::repeat_n(value, copies));
        }
        out
    }

    /// Computes `y = L_a x` for a function `x` on the vertices of the cut mesh, where
    ///
    /// ```text
    /// L_a = Σ_t ā_t · (cotan Laplacian of the single triangle t),   ā_t = (a_i + a_j + a_k)/3
    /// ```
    ///
    /// is the stiffness matrix of `∫ a ∇f·∇g` for piecewise linear `a`, `f`, `g` (SPEC §3.6,
    /// D11). `weight` holds `a` per vertex of the cut mesh; `None` means `a ≡ 1`, the cotan
    /// Laplacian `L` of the cut mesh (positive semidefinite). The two sides of the cut are not
    /// coupled.
    ///
    /// # Errors
    /// [`CutError::WrongLength`] if `x`, `y` or `weight` do not have
    /// [`CutTriangles::n_vertices`] entries.
    ///
    /// [formula: cut-weighted-laplacian]
    pub fn apply_laplacian(
        &self,
        weight: Option<&[f64]>,
        x: &[f64],
        y: &mut [f64],
    ) -> Result<(), CutError> {
        let n = self.n_vertices();
        check_length("x", "cut-mesh vertices", x.len(), n)?;
        check_length("y", "cut-mesh vertices", y.len(), n)?;
        if let Some(a) = weight {
            check_length("weight", "cut-mesh vertices", a.len(), n)?;
        }
        y.fill(0.0);
        for (face, cot) in self.faces.iter().zip(&self.cot) {
            let [i0, i1, i2] = face.map(|v| v as usize);
            let scale = match weight {
                Some(a) => 0.5 * (a[i0] + a[i1] + a[i2]) / 3.0,
                None => 0.5,
            };
            // The edge opposite corner k joins corners k + 1 and k + 2.
            for (k, (i, j)) in [(i1, i2), (i2, i0), (i0, i1)].into_iter().enumerate() {
                let flow = scale * cot[k] * (x[i] - x[j]);
                y[i] += flow;
                y[j] -= flow;
            }
        }
        Ok(())
    }

    /// The cut mesh as a validated [`TriMesh`] with a [`Geometry`] whose lengths are in the
    /// caller's unit (`scale` converts from the unit-area surface), and the maps of SPEC §7.
    ///
    /// # Errors
    /// [`CutError::DegenerateCutTriangle`] if the edge lengths of a triangle do not satisfy the
    /// strict triangle inequality in floating point arithmetic (a crossing extremely close to a
    /// vertex; the evaluation itself works with areas and inner products and is less strict).
    pub fn to_mesh(&self, scale: &UnitAreaScale) -> Result<CutMesh, CutError> {
        let mesh = TriMesh::new(self.n_vertices(), self.faces.clone())?;
        let mut lengths = vec![0.0f64; mesh.n_edges()];
        let mut seen = vec![false; mesh.n_edges()];
        for (t, len2) in self.len2.iter().enumerate() {
            for (k, squared) in len2.iter().enumerate() {
                let e = mesh.edge(3 * t + k);
                if !seen[e] {
                    seen[e] = true;
                    lengths[e] = scale.length_from_unit_area(squared.sqrt());
                }
            }
        }
        let geometry = Geometry::from_edge_lengths(&mesh, lengths).map_err(|e| match e {
            GeometryError::DegenerateFaces { faces } if !faces.is_empty() => {
                let first = &faces[0];
                CutError::DegenerateCutTriangle {
                    face: self.parent[first.face as usize],
                    vertices: first.vertices,
                    detail: format!(
                        "cut-mesh triangle {} with these cut-mesh vertices has the edge lengths \
                         {:?} ({} degenerate cut-mesh triangles in total)",
                        first.face,
                        first.lengths,
                        faces.len()
                    ),
                }
            }
            other => CutError::Geometry(other),
        })?;
        Ok(CutMesh {
            mesh,
            geometry,
            face_parent: self.parent.clone(),
            face_region: self.region.clone(),
            vertex_source: (0..self.n_vertices())
                .map(|v| self.vertex_source(v))
                .collect(),
            vertex_region: self.vertex_region.clone(),
        })
    }

    /// The halfedge of the input mesh that the edge `i -> j` of triangle `t` is, if `i` and `j`
    /// are both vertices of the input mesh.
    #[inline]
    pub(crate) fn original_halfedge(
        &self,
        mesh: &TriMesh,
        t: usize,
        i: u32,
        j: u32,
    ) -> Option<usize> {
        let f = self.parent[t] as usize;
        let face = mesh.faces()[f];
        (0..3)
            .find(|&k| face[k] == i && face[(k + 1) % 3] == j)
            .map(|k| 3 * f + k)
    }
}

impl CutMesh {
    /// Positions of the cut-mesh vertices in space, for a caller that has positions of the input
    /// vertices: an original vertex keeps its position, both copies of a cut point get
    /// `(1 − t)·x_a + t·x_b`, the three copies of a triple point its barycentric combination of
    /// the corners of its face.
    ///
    /// # Errors
    /// [`CutError::WrongLength`] if `positions` does not have one entry per vertex of `mesh`, or
    /// `cut` and `mesh` are not the ones this cut mesh was built from.
    pub fn positions(
        &self,
        cut: &Cut,
        mesh: &TriMesh,
        positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, CutError> {
        let points = cut.point_positions(mesh, positions)?;
        let triples = cut.triple_positions(mesh, positions)?;
        check_length(
            "cut points and input vertices",
            "cut-mesh vertices",
            positions.len() + 2 * points.len() + 3 * triples.len(),
            self.vertex_source.len(),
        )?;
        Ok(self
            .vertex_source
            .iter()
            .map(|source| match *source {
                CutVertex::Original(v) => positions[v as usize],
                CutVertex::Point(k) if (k as usize) < points.len() => points[k as usize],
                CutVertex::Point(k) => triples[k as usize - points.len()],
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dist2(p: [f64; 3], q: [f64; 3]) -> f64 {
        (0..3).map(|k| (p[k] - q[k]).powi(2)).sum()
    }

    fn point(x: [f64; 3], corners: [[f64; 3]; 3]) -> [f64; 3] {
        [0, 1, 2].map(|k| x[0] * corners[0][k] + x[1] * corners[1][k] + x[2] * corners[2][k])
    }

    #[test]
    fn barycentric_formulas_match_an_embedding() {
        let corners = [[0.1, 0.2, -0.3], [1.3, 0.4, 0.2], [0.5, 1.7, 0.9]];
        let len2 = [
            dist2(corners[0], corners[1]),
            dist2(corners[1], corners[2]),
            dist2(corners[2], corners[0]),
        ];
        let cross = {
            let u = [0, 1, 2].map(|k| corners[1][k] - corners[0][k]);
            let v = [0, 1, 2].map(|k| corners[2][k] - corners[0][k]);
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let face_area = 0.5 * dist2(cross, [0.0; 3]).sqrt();
        let x = [[0.2, 0.5, 0.3], [0.05, 0.15, 0.8], [0.7, 0.1, 0.2]];
        let embedded = x.map(|b| point(b, corners));
        // Squared lengths and inner products.
        for k in 0..3 {
            let (p, q) = (x[k], x[(k + 1) % 3]);
            let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
            let expect = dist2(embedded[k], embedded[(k + 1) % 3]);
            assert!((barycentric_dot(d, d, len2) - expect).abs() < 1e-14);
        }
        // Area and cotangents against the embedded triangle.
        let (area, sub_len2, cot) = subtriangle(x, len2, face_area).unwrap();
        let (a, b, c) = (sub_len2[0].sqrt(), sub_len2[1].sqrt(), sub_len2[2].sqrt());
        let heron = ddg_core::geometry::triangle_area(a, b, c).unwrap();
        assert!((area - heron).abs() < 1e-14, "{area} vs {heron}");
        for k in 0..3 {
            // The side opposite corner k is edge k + 1.
            let expect = ddg_core::geometry::corner_cotan(
                sub_len2[(k + 1) % 3].sqrt(),
                sub_len2[k].sqrt(),
                sub_len2[(k + 2) % 3].sqrt(),
            )
            .unwrap();
            assert!((cot[k] - expect).abs() < 1e-13, "{} vs {expect}", cot[k]);
        }
        // A reversed triangle has negative area and is rejected.
        assert!(subtriangle([x[0], x[2], x[1]], len2, face_area).is_none());
        // The whole face.
        let whole = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(subtriangle_area(whole, face_area), face_area);
        assert!((aspect_ratio([1.0, 1.0, 1.0], 3f64.sqrt() / 4.0) - 1.0).abs() < 1e-15);
    }
}
