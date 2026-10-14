//! The way back: a cut found on an intrinsic triangulation, expressed on the caller's polygons.
//!
//! The cut flow works on the export of the intrinsic triangulation: cut points on its edges,
//! triple points in its faces, straight segments across its faces. `ddg-intrinsic` carries them
//! back to the seam-cut triangle mesh: every cut point and every end of a segment becomes a point
//! of that mesh (inside a face, or on an edge), and every segment becomes a chain of straight
//! pieces that each lie in one face, joined at the points where the chain crosses an edge of the
//! mesh. A face of the seam-cut mesh is a triangle of the caller's triangulation, an edge is a
//! polygon edge or a diagonal, so those points and pieces are what [`CutOnPolygons`] describes:
//!
//! * a point on a polygon edge is an [`EdgePoint`] (an edge can carry several, ordered by `t`);
//! * a point on a diagonal is an [`InteriorPoint`] (a diagonal can carry several);
//! * a point inside a triangle is a [`JunctionPoint`] of any degree: a bend (two pieces), a
//!   triple point (three), and so on;
//! * a *path* runs between two points that are not crossings of a diagonal, through the
//!   crossings of diagonals: its consecutive points can lie in the same triangle, a triangle can
//!   hold several paths, and a path can enter and leave a triangle through the same side
//!   (through a bend inside it).
//!
//! **The rule for the points.** Every point keeps at least [`IntrinsicOptions::min_offset`]
//! (default `1e-6`) from the vertices and sides of its triangle: a point of an edge is moved
//! into `[min_offset, 1 − min_offset]`, a point inside a triangle has every barycentric
//! coordinate raised to `min_offset` and is renormalized. Moved points are counted in
//! [`CutOnPolygons::n_clamped_crossings`] (together with those the flow moved). A point at a
//! vertex of the mesh exactly is an error. The regions of the parts of a cut polygon do not come
//! from the regions of its vertices (a triangle with its three vertices in one region can be
//! crossed): they come from the sides of the oriented pieces, lower region id on the left.

use std::collections::HashMap;

use ddg_core::topology::NONE;
use ddg_intrinsic::{Segment, SurfacePoint};
use ddg_variational_cuts::Cut;

use super::describe::{
    CutLine, CutOnPolygons, EdgePoint, InteriorPoint, JunctionPoint, PathEnd, PolygonPath,
};
use super::error::{CuttingError, internal};
use super::session::{PreviewPolyline, Session};

/// A point of the seam-cut mesh.
#[derive(Debug, Clone, Copy)]
enum Place {
    Edge { edge: u32, t: f64 },
    Face { face: u32, bary: [f64; 3] },
}

/// A straight piece of the cut inside a face of the seam-cut mesh.
#[derive(Debug, Clone, Copy)]
struct Piece {
    a: u32,
    b: u32,
    face: u32,
    /// The regions on the left and on the right of `a -> b`.
    left: u32,
    right: u32,
}

impl Session {
    /// The cut on the caller's polygons for a session on an intrinsic triangulation.
    pub(super) fn describe_intrinsic(&self, cut: &Cut) -> Result<CutOnPolygons, CuttingError> {
        let t = self.tri.as_ref().expect("an intrinsic session");
        let m = &self.maps;
        let eps = t.options.min_offset;
        let n_points = cut.points.len();
        let export = &t.export;

        // Nodes of the cut on the triangulation, and its segments.
        let mut nodes: Vec<SurfacePoint> = cut
            .points
            .iter()
            .map(|p| export.edge_point(p.edge as usize, p.t))
            .collect();
        for triple in &cut.triple_points {
            nodes.push(SurfacePoint::Face {
                face: triple.face,
                bary: triple.barycentric,
            });
        }
        let segments: Vec<Segment> = cut
            .segments
            .iter()
            .map(|s| Segment {
                face: s.face,
                nodes: s.points,
            })
            .collect();
        let back = t
            .triangulation
            .segments_on_input(&nodes, &segments, &t.options.trace)?;

        // The regions on the two sides of every segment, in its direction.
        let sides = |s: usize| -> Result<(u32, u32), CuttingError> {
            let [p0, p1] = cut.segments[s].points;
            if p1 as usize >= n_points {
                let triple = &cut.triple_points[p1 as usize - n_points];
                let k = triple.points.iter().position(|&q| q == p0).ok_or_else(|| {
                    internal("a segment ends at a triple point that does not list its first end")
                })?;
                let face = export.mesh.faces()[triple.face as usize];
                Ok((
                    cut.vertex_region[face[k] as usize],
                    cut.vertex_region[face[(k + 1) % 3] as usize],
                ))
            } else {
                let [lo, hi] = cut.segment_regions(&export.mesh, s);
                Ok((lo, hi))
            }
        };

        // Places of the nodes and of the crossings, moved off vertices and sides.
        let mut n_clamped = cut.n_clamped;
        let mut places: Vec<Place> = Vec::with_capacity(back.nodes.len());
        let mut clamp = |point: SurfacePoint, what: &str| -> Result<Place, CuttingError> {
            match point {
                SurfacePoint::Vertex(v) => Err(CuttingError::Unsupported {
                    what: format!(
                        "{what} of the cut lies exactly at vertex {} of the seam-cut mesh (vertex \
                         {} of the polygon mesh)",
                        v, m.vertex_input[v as usize]
                    ),
                }),
                SurfacePoint::Edge { edge, t } => {
                    let c = t.clamp(eps, 1.0 - eps);
                    n_clamped += usize::from(c != t);
                    Ok(Place::Edge { edge, t: c })
                }
                SurfacePoint::Face { face, bary } => {
                    let raised = bary.map(|b| b.max(eps));
                    let sum: f64 = raised.iter().sum();
                    let c = raised.map(|b| b / sum);
                    n_clamped += usize::from(raised != bary);
                    Ok(Place::Face { face, bary: c })
                }
            }
        };
        for (i, &p) in back.nodes.iter().enumerate() {
            places.push(clamp(p, &format!("node {i}"))?);
        }
        let mut pieces: Vec<Piece> = Vec::new();
        for (s, path) in back.paths.iter().enumerate() {
            let (left, right) = sides(s)?;
            let [first, last] = cut.segments[s].points;
            let mut previous = first;
            for j in 0..path.faces.len() {
                let next = if j + 1 == path.faces.len() {
                    last
                } else {
                    places.push(clamp(path.points[j + 1], "a crossing of an edge")?);
                    (places.len() - 1) as u32
                };
                pieces.push(Piece {
                    a: previous,
                    b: next,
                    face: path.faces[j],
                    left,
                    right,
                });
                previous = next;
            }
        }

        // Points of the description.
        let mut description = CutOnPolygons {
            n_clamped_crossings: n_clamped,
            ..CutOnPolygons::default()
        };
        let x = &self.flow_positions;
        let mut edge_point_of = vec![NONE; places.len()];
        let mut junction_of = vec![NONE; places.len()];
        let mut interior_of: Vec<Option<InteriorPoint>> = vec![None; places.len()];
        for (g, place) in places.iter().enumerate() {
            match *place {
                Place::Edge { edge, t } => {
                    let e = edge as usize;
                    let [a, b] = self.flow_mesh.edge_vertices(e);
                    let position = [0, 1, 2].map(|k| (1.0 - t) * x[a][k] + t * x[b][k]);
                    if m.edge_is_polygon[e] {
                        let (va, vb) = (m.vertex_input[a], m.vertex_input[b]);
                        let tt = if va < vb { t } else { 1.0 - t };
                        let edge = [va.min(vb), va.max(vb)];
                        let side = m.used_seams.binary_search(&edge).ok().map(|_| {
                            let face = self.flow_mesh.face(self.flow_mesh.edge_halfedge(e));
                            m.triangle_polygon[m.face_triangle[face] as usize]
                        });
                        edge_point_of[g] = description.edge_points.len() as u32;
                        description.edge_points.push(EdgePoint {
                            edge,
                            t: tt,
                            position,
                            side,
                        });
                    } else {
                        let (c, d) = m.edge_corners[e];
                        let (lo, hi) = (c.min(d), c.max(d));
                        let tt = if m.corner_vertex[lo as usize] as usize == a {
                            t
                        } else {
                            1.0 - t
                        };
                        interior_of[g] = Some(InteriorPoint {
                            corners: [lo, hi],
                            t: tt,
                            position,
                        });
                    }
                }
                Place::Face { face, bary } => {
                    let f = self.flow_mesh.faces()[face as usize];
                    let position =
                        [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * x[f[c] as usize][k]).sum());
                    let triangle = m.face_triangle[face as usize];
                    junction_of[g] = description.junctions.len() as u32;
                    description.junctions.push(JunctionPoint {
                        position,
                        polygon: m.triangle_polygon[triangle as usize],
                        triangle,
                        corners: self.mesh.triangles[triangle as usize],
                        barycentric: bary,
                    });
                }
            }
        }

        // Pieces per polygon and node.
        let polygon_of =
            |piece: &Piece| m.triangle_polygon[m.face_triangle[piece.face as usize] as usize];
        let mut incident: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (i, piece) in pieces.iter().enumerate() {
            for g in [piece.a, piece.b] {
                incident.entry((polygon_of(piece), g)).or_default().push(i);
            }
        }
        for (g, place) in places.iter().enumerate() {
            if let Place::Face { face, .. } = place {
                // An inner point: every piece at it lies in its face.
                for p in pieces.iter().filter(|p| p.a == g as u32 || p.b == g as u32) {
                    if p.face != *face {
                        return Err(internal(format!(
                            "a piece at the point inside face {face} lies in face {}",
                            p.face
                        )));
                    }
                }
            }
        }
        let is_terminal = |g: u32| interior_of[g as usize].is_none();
        let mut order: Vec<usize> = (0..pieces.len()).collect();
        order.sort_by_key(|&i| (polygon_of(&pieces[i]), i));
        let mut used = vec![false; pieces.len()];
        let end_of = |g: u32| {
            if junction_of[g as usize] != NONE {
                PathEnd::Junction(junction_of[g as usize])
            } else {
                PathEnd::EdgePoint(edge_point_of[g as usize])
            }
        };
        for &s in &order {
            if used[s] {
                continue;
            }
            let polygon = polygon_of(&pieces[s]);
            let (a, b) = (pieces[s].a, pieces[s].b);
            let (start, forward) = if is_terminal(a) {
                (a, true)
            } else if is_terminal(b) {
                (b, false)
            } else {
                continue;
            };
            let mut walk = vec![start];
            let (mut current, mut from) = (s, start);
            loop {
                used[current] = true;
                let to = if pieces[current].a == from {
                    pieces[current].b
                } else {
                    pieces[current].a
                };
                walk.push(to);
                if is_terminal(to) {
                    break;
                }
                let list = &incident[&(polygon, to)];
                if list.len() != 2 {
                    return Err(CuttingError::Unsupported {
                        what: format!(
                            "a point on a diagonal in polygon {polygon} has {} pieces of the cut \
                             (a cut that runs along a diagonal, or a junction on one)",
                            list.len()
                        ),
                    });
                }
                current = if list[0] == current { list[1] } else { list[0] };
                from = to;
            }
            let (mut left, mut right) = if forward {
                (pieces[s].left, pieces[s].right)
            } else {
                (pieces[s].right, pieces[s].left)
            };
            if left > right {
                walk.reverse();
                std::mem::swap(&mut left, &mut right);
            }
            let interior = walk[1..walk.len() - 1]
                .iter()
                .map(|&g| interior_of[g as usize].clone().expect("a diagonal point"))
                .collect();
            description.paths.push(PolygonPath {
                polygon,
                start: end_of(walk[0]),
                interior,
                end: end_of(walk[walk.len() - 1]),
                regions: [left, right],
            });
        }
        if let Some(s) = used.iter().position(|u| !u) {
            return Err(CuttingError::Unsupported {
                what: format!(
                    "a piece of the cut in polygon {} lies along a diagonal of the triangulation",
                    polygon_of(&pieces[s])
                ),
            });
        }
        // Degrees: a point on a polygon edge has one piece in each polygon that uses the edge.
        description.paths.sort_by_key(|p| {
            let key = |end: PathEnd| match end {
                PathEnd::EdgePoint(i) => (0, i),
                PathEnd::Junction(i) => (1, i),
            };
            (p.polygon, key(p.start), key(p.end))
        });
        Ok(description)
    }

    /// The cut as lines for display, from the description ([`CutOnPolygons::lines`]).
    pub(super) fn preview_from_description(
        &self,
        cut: &Cut,
    ) -> Result<Vec<PreviewPolyline>, CuttingError> {
        let description = self.describe_intrinsic(cut)?;
        Ok(description
            .lines()
            .into_iter()
            .map(|l: CutLine| PreviewPolyline {
                positions: l.positions,
                closed: l.closed,
                regions: l.regions,
            })
            .collect())
    }
}
