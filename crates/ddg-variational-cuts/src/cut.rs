//! The cut: the zero set of the level-set function as points on mesh edges and segments across
//! faces (SPEC §3.2, §7).

use ddg_core::TriMesh;
use ddg_core::topology::NONE;

use crate::error::CutError;
use crate::problem::check_length;

/// A point of the cut on an edge of the input mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutPoint {
    /// The edge `[a, b]` (`a < b`) it lies on, as an edge id of the input mesh.
    pub edge: u32,
    /// Its parameter from `a`: the point is `(1 − t)·a + t·b`, with `0 < t < 1`.
    pub t: f64,
}

/// A straight piece of the cut across one face of the input mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CutSegment {
    /// Its two ends, as *nodes* of the cut: a node below `Cut::points.len()` is that cut point,
    /// a node `Cut::points.len() + m` is the triple point `m` ([`Cut::triple_points`]).
    ///
    /// In a face with two regions both ends are cut points; walking from the first to the
    /// second with the surface normal pointing up, the region with the lower id is on the left.
    /// In a face with three regions (SPEC §8.2) the first end is a cut point and the second the
    /// triple point of the face. The two regions a segment separates are those of its first
    /// end ([`Cut::segment_regions`]).
    pub points: [u32; 2],
    /// The face of the input mesh it crosses.
    pub face: u32,
}

/// A point where three regions meet, inside a face of the input mesh whose three vertices lie
/// in three regions (SPEC §8.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriplePoint {
    /// The face.
    pub face: u32,
    /// The cut points on its three edges: `points[k]` lies on the edge from corner `k` to
    /// corner `k + 1` of the face.
    pub points: [u32; 3],
    /// Its barycentric coordinates with respect to the three corners of the face: the centroid
    /// of the three cut points.
    ///
    /// [formula: cut-triple-point]
    pub barycentric: [f64; 3],
}

/// Where a curve of a cut network ends ([`CutCurve`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveEnd {
    /// On the boundary of the mesh.
    Boundary,
    /// At this triple point (an index into [`Cut::triple_points`]).
    Triple(u32),
    /// Nowhere: the curve is closed, its last point is joined to its first.
    Closed,
}

/// A curve of a cut network between two junctions, as returned by [`Cut::curves`]: a maximal
/// chain of cut points joined by segments, ending at triple points, on the boundary of the mesh,
/// or closed.
#[derive(Debug, Clone, PartialEq)]
pub struct CutCurve {
    /// The cut points of the curve, in order (indices into [`Cut::points`]).
    pub points: Vec<u32>,
    /// The polyline in space: the positions of the cut points, preceded by the position of the
    /// triple point the curve starts at and followed by that of the triple point it ends at,
    /// if any.
    pub positions: Vec<[f64; 3]>,
    /// The two regions it separates, `[lower id, higher id]`.
    pub regions: [u32; 2],
    /// Where it starts and where it ends.
    pub ends: [CurveEnd; 2],
}

/// One connected piece of a cut, as returned by [`Cut::polylines`].
#[derive(Debug, Clone, PartialEq)]
pub struct CutPolyline {
    /// The cut points of the piece, in order (indices into [`Cut::points`]).
    pub points: Vec<u32>,
    /// Their positions in space.
    pub positions: Vec<[f64; 3]>,
    /// Whether the piece is closed (the last point is joined to the first). An open piece ends
    /// on the boundary of the mesh at both ends.
    pub closed: bool,
}

/// The cut `γ` of a level-set function: closed polylines, or polylines ending on the boundary of
/// the mesh.
///
/// Points are numbered in order of their edge id, segments in order of their face id, so the
/// numbering depends only on the mesh and on which edges are crossed.
#[derive(Debug, Clone, PartialEq)]
pub struct Cut {
    /// The region of every vertex of the input mesh.
    pub vertex_region: Vec<u32>,
    /// The points of the cut.
    pub points: Vec<CutPoint>,
    /// The segments of the cut. A face with two regions has one segment, a face with three
    /// regions three (from its three cut points to its triple point, in the order of its edges).
    pub segments: Vec<CutSegment>,
    /// The triple points, in order of their face id. Empty for two regions.
    pub triple_points: Vec<TriplePoint>,
    /// For every edge of the input mesh, the index of its cut point, or
    /// [`NONE`](ddg_core::topology::NONE) if the edge is not crossed.
    pub edge_point: Vec<u32>,
    /// Number of crossings that were moved into `[δ, 1 − δ]` (SPEC §3.2, D9b).
    pub n_clamped: usize,
    /// For every cut point, whether its crossing was moved into `[δ, 1 − δ]`: the level-set
    /// function does not vanish where the cut point is.
    pub moved: Vec<bool>,
    /// The smallest `min(t, 1 − t)` over all crossings before they were moved; `0.5` if there is
    /// no crossing.
    pub smallest_offset: f64,
}

/// The parameter at which the linear interpolant of `φ` along an edge `[a, b]` with
/// `φ_a ≤ 0 < φ_b` or `φ_b ≤ 0 < φ_a` vanishes, measured from `a` (SPEC §3.2):
///
/// ```text
/// t = φ_a / (φ_a − φ_b)
/// ```
///
/// The result lies in `[0, 1]` for finite values of opposite sign class; it is `0` when
/// `φ_a = 0`.
///
/// [formula: cut-edge-crossing]
#[inline]
pub fn edge_crossing(phi_a: f64, phi_b: f64) -> f64 {
    let difference = phi_a - phi_b;
    if difference.is_finite() {
        phi_a / difference
    } else {
        // The difference of two huge values of opposite sign overflowed; halving is exact.
        (0.5 * phi_a) / (0.5 * phi_a - 0.5 * phi_b)
    }
}

/// The region of a vertex with level-set value `phi`: 0 where `φ ≤ 0`, 1 where `φ > 0`.
#[inline]
pub fn region_of(phi: f64) -> u32 {
    u32::from(phi > 0.0)
}

/// The region of every vertex for `N ≥ 2` level-set functions (SPEC §8.1): the index of the
/// smallest function at the vertex; ties go to the lowest index.
///
/// # Errors
/// [`CutError::UnsupportedRegionCount`] for fewer than two functions,
/// [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`].
pub fn regions_of(mesh: &TriMesh, phis: &[Vec<f64>]) -> Result<Vec<u32>, CutError> {
    if phis.len() < 2 {
        return Err(CutError::UnsupportedRegionCount { got: phis.len() });
    }
    for phi in phis {
        check_length(
            "level-set function",
            "vertices",
            phi.len(),
            mesh.n_vertices(),
        )?;
        if let Some(v) = phi.iter().position(|x| !x.is_finite()) {
            return Err(CutError::NonFiniteLevelSet {
                vertex: v as u32,
                value: phi[v],
            });
        }
    }
    Ok((0..mesh.n_vertices())
        .map(|v| {
            let mut best = 0;
            for k in 1..phis.len() {
                if phis[k][v] < phis[best][v] {
                    best = k;
                }
            }
            best as u32
        })
        .collect())
}

/// The corner of a face that is alone in its region, given the regions of its three corners:
/// `Ok(None)` if all three agree, `Ok(Some(k))` if corner `k` differs from the other two (which
/// agree), `Err(())` if all three differ.
#[inline]
pub(crate) fn lone_corner(regions: [u32; 3]) -> Result<Option<usize>, ()> {
    let [r0, r1, r2] = regions;
    if r0 == r1 && r1 == r2 {
        Ok(None)
    } else if r1 == r2 {
        Ok(Some(0))
    } else if r0 == r2 {
        Ok(Some(1))
    } else if r0 == r1 {
        Ok(Some(2))
    } else {
        Err(())
    }
}

impl Cut {
    /// Extracts the cut of the two-region level-set function `phi` (one value per vertex):
    /// region 0 is where `φ ≤ 0`, region 1 where `φ > 0`. A crossing closer than
    /// `min_crossing_offset` (a fraction of the edge) to a vertex is moved to that distance and
    /// counted in [`Cut::n_clamped`].
    ///
    /// Only signs and ratios of `phi` matter, so its length unit is irrelevant here.
    ///
    /// # Errors
    /// [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`], and
    /// [`CutError::InvalidOptions`] unless `0 ≤ min_crossing_offset < ½`.
    pub fn from_level_set(
        mesh: &TriMesh,
        phi: &[f64],
        min_crossing_offset: f64,
    ) -> Result<Cut, CutError> {
        if !(0.0..0.5).contains(&min_crossing_offset) {
            return Err(CutError::InvalidOptions {
                reason: format!(
                    "min_crossing_offset must lie in [0, 0.5), got {min_crossing_offset}"
                ),
            });
        }
        check_length(
            "level-set function",
            "vertices",
            phi.len(),
            mesh.n_vertices(),
        )?;
        if let Some(v) = phi.iter().position(|x| !x.is_finite()) {
            return Err(CutError::NonFiniteLevelSet {
                vertex: v as u32,
                value: phi[v],
            });
        }
        let vertex_region: Vec<u32> = phi.iter().map(|&x| region_of(x)).collect();
        Ok(Cut::assemble(
            mesh,
            vertex_region,
            |a, b| edge_crossing(phi[a], phi[b]),
            min_crossing_offset,
        ))
    }

    /// Extracts the cut network of `N ≥ 2` level-set functions (SPEC §8.1, §8.2): `phis[k]` has
    /// one value per vertex and is negative inside region `k`; the region of a vertex is the
    /// function that is smallest there (ties: the lowest `k`). An edge whose ends lie in the
    /// regions `i ≠ j` carries one cut point, where the interpolants of `φ_i` and `φ_j` are
    /// equal (paper Eq. 31),
    ///
    /// ```text
    /// t = ψ_a / (ψ_a − ψ_b),     ψ = φ_i − φ_j,
    /// ```
    ///
    /// moved into `[δ, 1 − δ]` like any crossing. A face with three regions gets a triple point
    /// at the centroid of its three cut points and three segments to it.
    ///
    /// With `phis = [φ, −φ]` this is [`Cut::from_level_set`] of `φ`, bit for bit.
    ///
    /// # Errors
    /// [`CutError::UnsupportedRegionCount`] for fewer than two functions,
    /// [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`], and
    /// [`CutError::InvalidOptions`] unless `0 ≤ min_crossing_offset < ½`.
    ///
    /// [formula: cut-edge-crossing-regions]
    pub fn from_regions(
        mesh: &TriMesh,
        phis: &[Vec<f64>],
        min_crossing_offset: f64,
    ) -> Result<Cut, CutError> {
        if !(0.0..0.5).contains(&min_crossing_offset) {
            return Err(CutError::InvalidOptions {
                reason: format!(
                    "min_crossing_offset must lie in [0, 0.5), got {min_crossing_offset}"
                ),
            });
        }
        let vertex_region = regions_of(mesh, phis)?;
        Ok(Cut::assemble(
            mesh,
            vertex_region.clone(),
            |a, b| {
                let (i, j) = (vertex_region[a] as usize, vertex_region[b] as usize);
                edge_crossing(phis[i][a] - phis[j][a], phis[i][b] - phis[j][b])
            },
            min_crossing_offset,
        ))
    }

    /// The cut for given regions of the vertices; `raw(a, b)` is the crossing parameter of the
    /// crossed edge `[a, b]` before it is moved.
    fn assemble(
        mesh: &TriMesh,
        vertex_region: Vec<u32>,
        raw: impl Fn(usize, usize) -> f64,
        min_crossing_offset: f64,
    ) -> Cut {
        let (lo, hi) = (min_crossing_offset, 1.0 - min_crossing_offset);
        let mut points = Vec::new();
        let mut edge_point = vec![NONE; mesh.n_edges()];
        let mut n_clamped = 0;
        let mut moved = Vec::new();
        let mut smallest_offset = 0.5f64;
        for (e, &[a, b]) in mesh.edges().iter().enumerate() {
            if vertex_region[a as usize] == vertex_region[b as usize] {
                continue;
            }
            let raw = raw(a as usize, b as usize);
            smallest_offset = smallest_offset.min(raw.min(1.0 - raw));
            let t = raw.clamp(lo, hi);
            if t != raw {
                n_clamped += 1;
            }
            moved.push(t != raw);
            edge_point[e] = points.len() as u32;
            points.push(CutPoint { edge: e as u32, t });
        }

        let mut segments = Vec::new();
        let mut triple_points: Vec<TriplePoint> = Vec::new();
        for (f, face) in mesh.faces().iter().enumerate() {
            let regions = face.map(|v| vertex_region[v as usize]);
            let k = match lone_corner(regions) {
                Ok(None) => continue,
                Ok(Some(k)) => k,
                Err(()) => {
                    // Three regions: the triple point is the centroid of the three cut points.
                    let on_edge = [0, 1, 2].map(|k| edge_point[mesh.edge(3 * f + k)]);
                    let mut barycentric = [0.0; 3];
                    for k in 0..3 {
                        let e = mesh.edge(3 * f + k);
                        let t = points[on_edge[k] as usize].t;
                        let s = if mesh.edges()[e][0] == face[k] {
                            t
                        } else {
                            1.0 - t
                        };
                        barycentric[k] += (1.0 - s) / 3.0;
                        barycentric[(k + 1) % 3] += s / 3.0;
                    }
                    let node = (points.len() + triple_points.len()) as u32;
                    for point in on_edge {
                        segments.push(CutSegment {
                            points: [point, node],
                            face: f as u32,
                        });
                    }
                    triple_points.push(TriplePoint {
                        face: f as u32,
                        points: on_edge,
                        barycentric,
                    });
                    continue;
                }
            };
            // The lone vertex a = corner k; p lies on the halfedge a -> b, q on c -> a. The
            // triangle (a, p, q) is counter-clockwise, so a's region is to the left of p -> q.
            let p = edge_point[mesh.edge(3 * f + k)];
            let q = edge_point[mesh.edge(3 * f + (k + 2) % 3)];
            let points = if regions[k] < regions[(k + 1) % 3] {
                [p, q]
            } else {
                [q, p]
            };
            segments.push(CutSegment {
                points,
                face: f as u32,
            });
        }
        Cut {
            vertex_region,
            points,
            segments,
            triple_points,
            edge_point,
            n_clamped,
            moved,
            smallest_offset,
        }
    }

    /// The two regions that segment `segment` separates, as `[lower id, higher id]`.
    ///
    /// # Panics
    /// If `segment` is out of range (it indexes like a slice).
    pub fn segment_regions(&self, mesh: &TriMesh, segment: usize) -> [u32; 2] {
        self.point_regions(mesh, self.segments[segment].points[0] as usize)
    }

    /// The positions of the triple points in space, for a caller that has vertex positions.
    ///
    /// # Errors
    /// As [`Cut::point_positions`].
    pub fn triple_positions(
        &self,
        mesh: &TriMesh,
        positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, CutError> {
        check_length("positions", "vertices", positions.len(), mesh.n_vertices())?;
        check_length(
            "the cut's edge map",
            "edges",
            self.edge_point.len(),
            mesh.n_edges(),
        )?;
        Ok(self
            .triple_points
            .iter()
            .map(|triple| {
                let face = mesh.faces()[triple.face as usize];
                [0, 1, 2].map(|k| {
                    (0..3)
                        .map(|c| triple.barycentric[c] * positions[face[c] as usize][k])
                        .sum()
                })
            })
            .collect())
    }

    /// The cut network as curves in space (SPEC §8.7): every maximal chain of cut points joined
    /// by segments between two cut points, with the two regions it separates and what it ends
    /// at. For two regions these are the pieces of [`Cut::polylines`].
    ///
    /// The order is fixed by the numbering of the cut points: first the curves that are not
    /// closed, in the order of their lower-numbered end point, each starting there; then the
    /// closed ones, in the order of their lowest-numbered point, starting there and going on to
    /// the lower-numbered of its two neighbours. A curve that consists of one cut point between
    /// two triple points starts at the triple point that comes first in the segments of the cut.
    ///
    /// # Errors
    /// As [`Cut::point_positions`].
    pub fn curves(
        &self,
        mesh: &TriMesh,
        positions: &[[f64; 3]],
    ) -> Result<Vec<CutCurve>, CutError> {
        let at = self.point_positions(mesh, positions)?;
        let triple_at = self.triple_positions(mesh, positions)?;
        let n = self.points.len();
        let mut neighbours = vec![[NONE; 2]; n];
        let mut triples = vec![[NONE; 2]; n];
        for segment in &self.segments {
            let [p, q] = segment.points;
            if (q as usize) < n {
                for (a, b) in [(p, q), (q, p)] {
                    let slot = usize::from(neighbours[a as usize][0] != NONE);
                    neighbours[a as usize][slot] = b;
                }
            } else {
                let slot = usize::from(triples[p as usize][0] != NONE);
                triples[p as usize][slot] = q - n as u32;
            }
        }
        let end_of = |triple: u32| {
            if triple == NONE {
                CurveEnd::Boundary
            } else {
                CurveEnd::Triple(triple)
            }
        };
        let mut visited = vec![false; n];
        let mut out = Vec::new();
        for closed in [false, true] {
            for start in 0..n {
                let [a, b] = neighbours[start];
                if visited[start] || (b == NONE) == closed {
                    continue;
                }
                let mut points = vec![start as u32];
                visited[start] = true;
                let (mut previous, mut current) = (start as u32, if closed { a.min(b) } else { a });
                while current != NONE && !visited[current as usize] {
                    visited[current as usize] = true;
                    points.push(current);
                    let [x, y] = neighbours[current as usize];
                    let next = if x == previous { y } else { x };
                    previous = current;
                    current = next;
                }
                let last = points[points.len() - 1] as usize;
                let ends = if closed {
                    [CurveEnd::Closed; 2]
                } else if points.len() == 1 {
                    [end_of(triples[start][0]), end_of(triples[start][1])]
                } else {
                    [end_of(triples[start][0]), end_of(triples[last][0])]
                };
                let mut line = Vec::with_capacity(points.len() + 2);
                if let CurveEnd::Triple(m) = ends[0] {
                    line.push(triple_at[m as usize]);
                }
                line.extend(points.iter().map(|&p| at[p as usize]));
                if let CurveEnd::Triple(m) = ends[1] {
                    line.push(triple_at[m as usize]);
                }
                out.push(CutCurve {
                    regions: self.point_regions(mesh, start),
                    points,
                    positions: line,
                    ends,
                });
            }
        }
        Ok(out)
    }

    /// The two regions that meet at cut point `point`, as `[lower id, higher id]`.
    ///
    /// # Panics
    /// If `point` is out of range (it indexes like a slice).
    pub fn point_regions(&self, mesh: &TriMesh, point: usize) -> [u32; 2] {
        let [a, b] = mesh.edge_vertices(self.points[point].edge as usize);
        let (ra, rb) = (self.vertex_region[a], self.vertex_region[b]);
        [ra.min(rb), ra.max(rb)]
    }

    /// For every vertex of the input mesh, whether it is an end of a crossed edge (the
    /// *cut-adjacent* vertices of SPEC §3.2).
    pub fn cut_adjacent_vertices(&self, mesh: &TriMesh) -> Vec<bool> {
        let mut adjacent = vec![false; self.vertex_region.len()];
        for point in &self.points {
            let [a, b] = mesh.edge_vertices(point.edge as usize);
            adjacent[a] = true;
            adjacent[b] = true;
        }
        adjacent
    }

    /// The positions of the cut points in space, for a caller that has vertex positions:
    /// `(1 − t)·x_a + t·x_b`.
    ///
    /// # Errors
    /// [`CutError::WrongLength`] if `positions` does not have one entry per vertex, or if the
    /// cut was not extracted from `mesh`.
    pub fn point_positions(
        &self,
        mesh: &TriMesh,
        positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, CutError> {
        check_length("positions", "vertices", positions.len(), mesh.n_vertices())?;
        check_length(
            "the cut's edge map",
            "edges",
            self.edge_point.len(),
            mesh.n_edges(),
        )?;
        Ok(self
            .points
            .iter()
            .map(|point| {
                let [a, b] = mesh.edge_vertices(point.edge as usize);
                let (xa, xb) = (positions[a], positions[b]);
                [0, 1, 2].map(|k| (1.0 - point.t) * xa[k] + point.t * xb[k])
            })
            .collect())
    }

    /// The cut as polylines in space: every connected piece of the cut with its points in
    /// order, and their positions `(1 − t)·x_a + t·x_b`.
    ///
    /// A piece is either closed or ends on the boundary of the mesh at both ends. The order is
    /// fixed by the numbering of the cut points: the open pieces come first, in the order of
    /// their lower-numbered end, and start there; then the closed pieces, in the order of their
    /// lowest-numbered point, starting there and going on to the lower-numbered of its two
    /// neighbours. A closed piece does not repeat its first point at the end.
    ///
    /// This is for two regions. On a network with triple points the pieces stop at the last cut
    /// point before a triple point; use [`Cut::curves`] there.
    ///
    /// # Errors
    /// As [`Cut::point_positions`].
    pub fn polylines(
        &self,
        mesh: &TriMesh,
        positions: &[[f64; 3]],
    ) -> Result<Vec<CutPolyline>, CutError> {
        let at = self.point_positions(mesh, positions)?;
        let n = self.points.len();
        // Every cut point lies on one or two segments.
        let mut neighbours = vec![[NONE; 2]; n];
        for segment in &self.segments {
            let [p, q] = segment.points;
            if q as usize >= n {
                // A segment to a triple point: see `Cut::curves`.
                continue;
            }
            for (a, b) in [(p, q), (q, p)] {
                let slot = usize::from(neighbours[a as usize][0] != NONE);
                neighbours[a as usize][slot] = b;
            }
        }
        let mut visited = vec![false; n];
        let mut out = Vec::new();
        let walk = |start: usize, first: u32, visited: &mut Vec<bool>| -> Vec<u32> {
            let mut points = vec![start as u32];
            visited[start] = true;
            let (mut previous, mut current) = (start as u32, first);
            while current != NONE && !visited[current as usize] {
                visited[current as usize] = true;
                points.push(current);
                let [a, b] = neighbours[current as usize];
                let next = if a == previous { b } else { a };
                previous = current;
                current = next;
            }
            points
        };
        for closed in [false, true] {
            for start in 0..n {
                let [a, b] = neighbours[start];
                if visited[start] || a == NONE || (b == NONE) == closed {
                    continue;
                }
                let first = if closed { a.min(b) } else { a };
                let points = walk(start, first, &mut visited);
                out.push(CutPolyline {
                    positions: points.iter().map(|&p| at[p as usize]).collect(),
                    points,
                    closed,
                });
            }
        }
        Ok(out)
    }

    /// Interpolates a per-vertex function of the input mesh linearly to the nodes of the cut:
    /// the cut points, then the triple points.
    pub(crate) fn interpolate_nodes(&self, mesh: &TriMesh, values: &[f64]) -> Vec<f64> {
        let mut out = self.interpolate(mesh, values);
        out.extend(self.triple_points.iter().map(|triple| {
            let face = mesh.faces()[triple.face as usize];
            (0..3)
                .map(|c| triple.barycentric[c] * values[face[c] as usize])
                .sum::<f64>()
        }));
        out
    }

    /// Interpolates a per-vertex function of the input mesh linearly to the cut points.
    pub(crate) fn interpolate(&self, mesh: &TriMesh, values: &[f64]) -> Vec<f64> {
        self.points
            .iter()
            .map(|point| {
                let [a, b] = mesh.edge_vertices(point.edge as usize);
                (1.0 - point.t) * values[a] + point.t * values[b]
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_parameter() {
        assert_eq!(edge_crossing(-1.0, 3.0), 0.25);
        assert_eq!(edge_crossing(3.0, -1.0), 0.75);
        assert_eq!(edge_crossing(0.0, 2.0), 0.0);
        assert_eq!(edge_crossing(-0.0, 2.0), 0.0);
        assert_eq!(edge_crossing(5.0, 0.0), 1.0);
        // Overflow of the difference is avoided.
        assert_eq!(edge_crossing(-f64::MAX, f64::MAX), 0.5);
        assert!((edge_crossing(f64::MAX, -0.5 * f64::MAX) - 2.0 / 3.0).abs() < 1e-15);
        // Tiny values keep their ratio.
        assert!((edge_crossing(-1e-300, 3e-300) - 0.25).abs() < 1e-15);
    }

    #[test]
    fn lone_corner_cases() {
        assert_eq!(lone_corner([1, 1, 1]), Ok(None));
        assert_eq!(lone_corner([0, 1, 1]), Ok(Some(0)));
        assert_eq!(lone_corner([1, 0, 1]), Ok(Some(1)));
        assert_eq!(lone_corner([1, 1, 0]), Ok(Some(2)));
        assert_eq!(lone_corner([0, 1, 2]), Err(()));
    }

    #[test]
    fn square_cut_through_two_faces() {
        // Unit square 0-1-2-3 split along 0-2; phi separates vertex 1 from the rest.
        let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
        let cut = Cut::from_level_set(&mesh, &[-1.0, 3.0, -1.0, -2.0], 1e-3).unwrap();
        assert_eq!(cut.vertex_region, vec![0, 1, 0, 0]);
        // Edges: [0,1], [0,2], [0,3], [1,2], [2,3]; crossed: [0,1] and [1,2].
        assert_eq!(cut.edge_point, vec![0, NONE, NONE, 1, NONE]);
        assert_eq!(
            cut.points,
            vec![CutPoint { edge: 0, t: 0.25 }, CutPoint { edge: 3, t: 0.75 }]
        );
        // Face 0 = (0, 1, 2): the lone vertex is 1 (region 1); p on 1 -> 2 (point 1), q on
        // 0 -> 1 (point 0). Region 1 is left of p -> q, so region 0 is left of q -> p.
        assert_eq!(
            cut.segments,
            vec![CutSegment {
                points: [0, 1],
                face: 0
            }]
        );
        assert_eq!(cut.n_clamped, 0);
        assert_eq!(cut.smallest_offset, 0.25);
        assert_eq!(cut.point_regions(&mesh, 0), [0, 1]);
        assert_eq!(
            cut.cut_adjacent_vertices(&mesh),
            vec![true, true, true, false]
        );
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let p = cut.point_positions(&mesh, &positions).unwrap();
        assert_eq!(p, vec![[0.25, 0.0, 0.0], [1.0, 0.75, 0.0]]);
    }

    #[test]
    fn zero_at_a_vertex_is_clamped_and_counted() {
        let mesh = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap();
        let cut = Cut::from_level_set(&mesh, &[0.0, 1.0, 0.0, -1.0], 0.01).unwrap();
        // Vertices 0 and 2 have phi = 0: region 0. Crossed: [0,1] at t = 0, [1,2] at t = 1.
        assert_eq!(cut.vertex_region, vec![0, 1, 0, 0]);
        assert_eq!(cut.points[0].t, 0.01);
        assert_eq!(cut.points[1].t, 0.99);
        assert_eq!(cut.n_clamped, 2);
        assert_eq!(cut.smallest_offset, 0.0);
    }
}
