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
    /// Its two end points, as indices into [`Cut::points`]. Walking from the first to the second
    /// with the surface normal pointing up, the region with the lower id is on the left.
    pub points: [u32; 2],
    /// The face of the input mesh it crosses.
    pub face: u32,
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
    /// The segments of the cut. A face has zero or one segment.
    pub segments: Vec<CutSegment>,
    /// For every edge of the input mesh, the index of its cut point, or
    /// [`NONE`](ddg_core::topology::NONE) if the edge is not crossed.
    pub edge_point: Vec<u32>,
    /// Number of crossings that were moved into `[δ, 1 − δ]` (SPEC §3.2, D9b).
    pub n_clamped: usize,
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
        if !(min_crossing_offset >= 0.0 && min_crossing_offset < 0.5) {
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

        let (lo, hi) = (min_crossing_offset, 1.0 - min_crossing_offset);
        let mut points = Vec::new();
        let mut edge_point = vec![NONE; mesh.n_edges()];
        let mut n_clamped = 0;
        let mut smallest_offset = 0.5f64;
        for (e, &[a, b]) in mesh.edges().iter().enumerate() {
            if vertex_region[a as usize] == vertex_region[b as usize] {
                continue;
            }
            let raw = edge_crossing(phi[a as usize], phi[b as usize]);
            smallest_offset = smallest_offset.min(raw.min(1.0 - raw));
            let t = raw.clamp(lo, hi);
            if t != raw {
                n_clamped += 1;
            }
            edge_point[e] = points.len() as u32;
            points.push(CutPoint { edge: e as u32, t });
        }

        let mut segments = Vec::new();
        for (f, face) in mesh.faces().iter().enumerate() {
            let regions = face.map(|v| vertex_region[v as usize]);
            let k = match lone_corner(regions) {
                Ok(None) => continue,
                Ok(Some(k)) => k,
                Err(()) => {
                    return Err(CutError::UnsupportedTriplePoint {
                        face: f as u32,
                        vertices: *face,
                        regions,
                    });
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
        Ok(Cut {
            vertex_region,
            points,
            segments,
            edge_point,
            n_clamped,
            smallest_offset,
        })
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
            vec![
                CutPoint { edge: 0, t: 0.25 },
                CutPoint { edge: 3, t: 0.75 }
            ]
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
