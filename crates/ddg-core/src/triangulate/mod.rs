//! Triangulation of polygons the way Blender does it.
//!
//! The algorithms of this workspace run on triangles, but meshes arrive as polygons and results
//! must be reported per original polygon corner. [`triangulate_polygons`] splits every polygon
//! exactly like Blender's own tessellation (the one behind `Mesh.loop_triangles`), so the
//! triangles used here are the triangles Blender draws:
//!
//! * a triangle is passed through;
//! * a quad is split along the diagonal from its first to its third corner, unless that would
//!   fold the two triangles over each other, in which case the other diagonal is used;
//! * a polygon with five or more corners is projected along its Newell normal and filled by
//!   Blender's 2D ear clipping.
//!
//! The caller's polygons are never modified. Every triangle corner records the polygon corner
//! ("loop") it came from, and every triangle side records whether it is an original polygon
//! edge or a diagonal introduced here — later stages must not place cuts or seams on diagonals.
//!
//! # Agreement with Blender
//!
//! `scripts/blender-tessellation-diff.sh` compares this module with a real Blender on random
//! quads and n-gons; see the crate README for the measured agreement. Decisions are made in
//! `f32` because Blender stores positions in single precision: positions are converted with
//! `as f32`, which is exact for values that came from Blender.
//!
//! # License
//!
//! The splitting rules live in the private `blender` submodule, a port of Blender source code
//! under the **GNU GPL, version 2 or later**. That submodule, and therefore this module, is
//! compiled only with the cargo feature `blender-tessellation` (enabled by default). Building
//! the crate without that feature leaves only MIT-licensed code.

mod blender;

use std::fmt;

use crate::topology::NONE;

/// Why a polygon soup could not be triangulated.
#[derive(Debug, Clone, PartialEq)]
pub enum TriangulateError {
    /// There are corners but no polygon starts, or more corners than `u32` ids can address.
    BadCounts {
        /// Number of corners.
        n_corners: usize,
        /// Number of polygons.
        n_faces: usize,
    },
    /// `face_starts` is not a valid list of polygon offsets.
    BadFaceStart {
        /// The polygon whose start is wrong.
        face: u32,
        /// Its start offset.
        start: u32,
        /// The start of the preceding polygon (0 for the first polygon).
        previous_start: u32,
        /// Total number of corners.
        n_corners: usize,
    },
    /// A polygon has fewer than three corners.
    FaceTooSmall {
        /// The polygon.
        face: u32,
        /// Its number of corners.
        size: u32,
    },
    /// A corner references a vertex that does not exist.
    VertexIndexOutOfRange {
        /// The polygon.
        face: u32,
        /// The corner (index into `corner_verts`).
        corner: u32,
        /// The offending vertex index.
        index: u32,
        /// Number of positions.
        n_vertices: usize,
    },
    /// A position is NaN, infinite, or too large for single precision.
    NonFinitePosition {
        /// The vertex.
        vertex: u32,
        /// Its position.
        position: [f64; 3],
    },
}

impl fmt::Display for TriangulateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TriangulateError::BadCounts { n_corners, n_faces } => {
                write!(
                    f,
                    "invalid polygon soup: {n_corners} corners in {n_faces} polygons"
                )
            }
            TriangulateError::BadFaceStart {
                face,
                start,
                previous_start,
                n_corners,
            } => {
                if *face == 0 {
                    write!(
                        f,
                        "polygon 0 starts at corner {start}, but the first polygon must start at corner 0"
                    )
                } else if (*start as usize) > *n_corners {
                    write!(
                        f,
                        "polygon {face} starts at corner {start}, but there are only {n_corners} corners"
                    )
                } else {
                    write!(
                        f,
                        "polygon {face} starts at corner {start}, before polygon {}, which starts at corner \
                         {previous_start}: polygon starts must not decrease",
                        face - 1
                    )
                }
            }
            TriangulateError::FaceTooSmall { face, size } => {
                write!(
                    f,
                    "polygon {face} has {size} corners, at least 3 are required"
                )
            }
            TriangulateError::VertexIndexOutOfRange {
                face,
                corner,
                index,
                n_vertices,
            } => write!(
                f,
                "corner {corner} (polygon {face}) references vertex {index}, but there are {n_vertices} positions"
            ),
            TriangulateError::NonFinitePosition {
                vertex,
                position: [x, y, z],
            } => write!(
                f,
                "vertex {vertex} has position ({x}, {y}, {z}), which is not finite in single precision"
            ),
        }
    }
}

impl std::error::Error for TriangulateError {}

/// The triangles of a polygon soup and their relation to the polygons.
///
/// Triangles are grouped by polygon, in polygon order: polygon `p` with `n` corners starting at
/// corner `s` owns the `n − 2` triangles starting at triangle `s − 2 p`. Triangles keep the
/// orientation of their polygon. All arrays have one entry per triangle.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Triangulation {
    /// The vertex indices of each triangle; ready for [`crate::topology::TriMesh::new`].
    pub triangles: Vec<[u32; 3]>,
    /// For each triangle corner, the polygon corner (Blender: loop) it came from, as an index
    /// into `corner_verts`. `triangles[t][k] == corner_verts[triangle_corners[t][k]]`.
    /// This is Blender's `loop_triangles[t].loops`.
    pub triangle_corners: Vec<[u32; 3]>,
    /// For each triangle, the polygon it belongs to.
    pub triangle_polygons: Vec<u32>,
    /// For each triangle side, whether it is a diagonal introduced by triangulation (`true`) or
    /// an edge of the original polygon (`false`). Side `k` runs from corner `k` to corner
    /// `(k + 1) % 3`, i.e. it is halfedge `3 t + k` of a mesh built from `triangles`.
    pub is_diagonal: Vec<[bool; 3]>,
}

impl Triangulation {
    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Copies per-triangle-corner values to the original polygon corners:
    /// `out[triangle_corners[t][k]] = values[t][k]`. `n_corners` is the length of the
    /// `corner_verts` array that was triangulated.
    ///
    /// A polygon corner shared by several triangles receives the value of the last of them; for
    /// data that is continuous across a polygon's diagonals (per-vertex data, or UVs when no
    /// seam runs along a diagonal) all of them agree. A corner that no triangle uses keeps
    /// `fill` (this happens only if `n_corners` exceeds the triangulated corner count).
    ///
    /// Returns `None` if `values` does not have one entry per triangle or if `n_corners` is
    /// smaller than the number of corners that were triangulated.
    pub fn to_polygon_corners<T: Copy>(
        &self,
        values: &[[T; 3]],
        n_corners: usize,
        fill: T,
    ) -> Option<Vec<T>> {
        if values.len() != self.triangles.len() {
            return None;
        }
        let mut out = vec![fill; n_corners];
        for (corners, value) in self.triangle_corners.iter().zip(values) {
            for k in 0..3 {
                *out.get_mut(corners[k] as usize)? = value[k];
            }
        }
        Some(out)
    }
}

/// Triangulates polygons given in Blender's layout, matching Blender's tessellation.
///
/// * `positions` — vertex positions;
/// * `corner_verts` — for every polygon corner (Blender: loop), its vertex index
///   (`MeshLoop.vertex_index`);
/// * `face_starts` — for every polygon, the index of its first corner
///   (`MeshPolygon.loop_start`). Starts must begin at 0 and increase by at least 3; the last
///   polygon extends to the end of `corner_verts`.
///
/// Vertices are not renumbered and `positions` may contain vertices that no polygon uses.
/// Polygons with repeated vertices are triangulated like any other; the resulting triangles
/// then fail [`crate::topology::TriMesh::new`], which names the triangle.
///
/// # Errors
/// See [`TriangulateError`]. Nothing is triangulated unless the whole input is valid.
pub fn triangulate_polygons(
    positions: &[[f64; 3]],
    corner_verts: &[u32],
    face_starts: &[u32],
) -> Result<Triangulation, TriangulateError> {
    let n_corners = corner_verts.len();
    let n_faces = face_starts.len();
    if n_corners >= NONE as usize || (n_faces == 0 && n_corners != 0) {
        return Err(TriangulateError::BadCounts { n_corners, n_faces });
    }
    // Polygon extents.
    let mut previous = 0u32;
    for (f, &start) in face_starts.iter().enumerate() {
        if (f == 0 && start != 0) || start < previous || start as usize > n_corners {
            return Err(TriangulateError::BadFaceStart {
                face: f as u32,
                start,
                previous_start: previous,
                n_corners,
            });
        }
        previous = start;
    }
    for (f, &start) in face_starts.iter().enumerate() {
        let end = face_starts.get(f + 1).map_or(n_corners, |&e| e as usize);
        let size = end - start as usize;
        if size < 3 {
            return Err(TriangulateError::FaceTooSmall {
                face: f as u32,
                size: size as u32,
            });
        }
    }
    for (f, &start) in face_starts.iter().enumerate() {
        let end = face_starts.get(f + 1).map_or(n_corners, |&e| e as usize);
        for (offset, &index) in corner_verts[start as usize..end].iter().enumerate() {
            if index as usize >= positions.len() {
                return Err(TriangulateError::VertexIndexOutOfRange {
                    face: f as u32,
                    corner: start + offset as u32,
                    index,
                    n_vertices: positions.len(),
                });
            }
        }
    }
    // Only positions that are used need to be representable.
    let single: Vec<[f32; 3]> = positions.iter().map(|p| p.map(|x| x as f32)).collect();
    for &v in corner_verts {
        if !single[v as usize].iter().all(|x| x.is_finite()) {
            return Err(TriangulateError::NonFinitePosition {
                vertex: v,
                position: positions[v as usize],
            });
        }
    }

    let n_triangles = n_corners - 2 * n_faces;
    let mut out = Triangulation {
        triangles: Vec::with_capacity(n_triangles),
        triangle_corners: Vec::with_capacity(n_triangles),
        triangle_polygons: Vec::with_capacity(n_triangles),
        is_diagonal: Vec::with_capacity(n_triangles),
    };
    let mut scratch = blender::PolyFillScratch::default();
    let mut corner_positions: Vec<[f32; 3]> = Vec::new();
    let mut local: Vec<[u32; 3]> = Vec::new();
    for (f, &start) in face_starts.iter().enumerate() {
        let end = face_starts.get(f + 1).map_or(n_corners, |&e| e as usize);
        let size = (end - start as usize) as u32;
        corner_positions.clear();
        corner_positions.extend(
            corner_verts[start as usize..end]
                .iter()
                .map(|&v| single[v as usize]),
        );
        local.clear();
        blender::tessellate_face(&corner_positions, &mut scratch, &mut local);
        debug_assert_eq!(local.len() as u32, size - 2);
        for tri in &local {
            let corners = tri.map(|i| start + i);
            out.triangles
                .push(corners.map(|c| corner_verts[c as usize]));
            out.triangle_corners.push(corners);
            out.triangle_polygons.push(f as u32);
            out.is_diagonal.push([0, 1, 2].map(|k| {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                // A polygon edge joins cyclically consecutive corners.
                !((a + 1) % size == b || (b + 1) % size == a)
            }));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tri_set(t: &Triangulation) -> Vec<[u32; 3]> {
        t.triangle_corners.clone()
    }

    #[test]
    fn triangles_pass_through() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ];
        let t = triangulate_polygons(&p, &[0, 1, 2, 2, 1, 3], &[0, 3]).unwrap();
        assert_eq!(t.triangles, vec![[0, 1, 2], [2, 1, 3]]);
        assert_eq!(t.triangle_corners, vec![[0, 1, 2], [3, 4, 5]]);
        assert_eq!(t.triangle_polygons, vec![0, 1]);
        assert_eq!(t.is_diagonal, vec![[false; 3]; 2]);
        assert_eq!(t.n_triangles(), 2);
    }

    #[test]
    fn convex_quad_uses_first_third_diagonal() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let t = triangulate_polygons(&p, &[0, 1, 2, 3], &[0]).unwrap();
        assert_eq!(tri_set(&t), vec![[0, 1, 2], [0, 2, 3]]);
        // Sides: (0,1) edge, (1,2) edge, (2,0) diagonal; (0,2) diagonal, (2,3) edge, (3,0) edge.
        assert_eq!(
            t.is_diagonal,
            vec![[false, false, true], [true, false, false]]
        );
    }

    #[test]
    fn concave_quad_flips_the_diagonal() {
        // A dart whose reflex corner is corner 2: the diagonal 0-2 lies inside, no flip.
        let p = [
            [0.0, 0.0, 0.0],
            [2.0, -1.0, 0.0],
            [0.5, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let t = triangulate_polygons(&p, &[0, 1, 2, 3], &[0]).unwrap();
        assert_eq!(tri_set(&t), vec![[0, 1, 2], [0, 2, 3]]);
        // The same dart starting one corner later: reflex corner is now corner 1, and the
        // diagonal 0-2 (vertices 1-3) lies outside the polygon. Blender flips to 1-3.
        let t = triangulate_polygons(&p, &[1, 2, 3, 0], &[0]).unwrap();
        assert_eq!(tri_set(&t), vec![[0, 1, 3], [1, 2, 3]]);
        assert_eq!(t.triangles, vec![[1, 2, 0], [2, 3, 0]]);
        assert_eq!(
            t.is_diagonal,
            vec![[false, true, false], [false, false, true]]
        );
    }

    #[test]
    fn convex_pentagon_is_clipped_evenly() {
        // Regular pentagon in the XY plane, counter-clockwise seen from +Z.
        let p: Vec<[f64; 3]> = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let t = triangulate_polygons(&p, &[0, 1, 2, 3, 4], &[0]).unwrap();
        assert_eq!(t.n_triangles(), 3);
        // Blender's even clipping: ear at corner 0, then at corner 2, then the rest.
        assert_eq!(tri_set(&t), vec![[4, 0, 1], [1, 2, 3], [1, 3, 4]]);
        assert_eq!(t.triangle_polygons, vec![0, 0, 0]);
        assert_eq!(
            t.is_diagonal,
            vec![
                [false, false, true],
                [false, false, true],
                [true, false, true]
            ]
        );
    }

    #[test]
    fn validation_errors() {
        let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let err =
            |cv: &[u32], fs: &[u32]| triangulate_polygons(&p, cv, fs).unwrap_err().to_string();
        assert_eq!(
            err(&[0, 1, 2], &[]),
            "invalid polygon soup: 3 corners in 0 polygons"
        );
        assert_eq!(
            err(&[0, 1, 2], &[1]),
            "polygon 0 starts at corner 1, but the first polygon must start at corner 0"
        );
        assert_eq!(
            err(&[0, 1], &[0]),
            "polygon 0 has 2 corners, at least 3 are required"
        );
        assert_eq!(
            err(&[0, 1, 2, 0, 1], &[0, 3]),
            "polygon 1 has 2 corners, at least 3 are required"
        );
        assert_eq!(
            err(&[0, 1, 2, 0, 1, 2], &[0, 2]),
            "polygon 0 has 2 corners, at least 3 are required"
        );
        assert_eq!(
            err(&[0, 1, 2, 0, 1, 2], &[0, 9]),
            "polygon 1 starts at corner 9, but there are only 6 corners"
        );
        assert_eq!(
            err(&[0, 1, 2, 0, 1, 2, 0, 1, 2], &[0, 6, 3]),
            "polygon 2 starts at corner 3, before polygon 1, which starts at corner 6: polygon starts must not \
             decrease"
        );
        assert_eq!(
            err(&[0, 1, 3], &[0]),
            "corner 2 (polygon 0) references vertex 3, but there are 3 positions"
        );
        let bad = [
            [0.0, 0.0, 0.0],
            [1.0, f64::NAN, 0.0],
            [0.0, 1.0, 0.0],
            [1e300, 0.0, 0.0],
        ];
        assert_eq!(
            triangulate_polygons(&bad, &[0, 1, 2], &[0])
                .unwrap_err()
                .to_string(),
            "vertex 1 has position (1, NaN, 0), which is not finite in single precision"
        );
        let e = triangulate_polygons(&bad, &[0, 3, 2], &[0]).unwrap_err();
        assert!(
            matches!(e, TriangulateError::NonFinitePosition { vertex: 3, .. }),
            "{e:?}"
        );
        assert!(
            e.to_string()
                .ends_with("which is not finite in single precision")
        );
        // An unused bad vertex does not matter.
        assert!(triangulate_polygons(&bad, &[0, 2, 0], &[0]).is_ok());
        // The empty soup is fine.
        assert_eq!(
            triangulate_polygons(&p, &[], &[]).unwrap(),
            Triangulation::default()
        );
    }

    #[test]
    fn to_polygon_corners_maps_back() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let t = triangulate_polygons(&p, &[0, 1, 2, 3], &[0]).unwrap();
        // Per-triangle-corner values that are really per-vertex: 10 * vertex.
        let values: Vec<[u32; 3]> = t.triangles.iter().map(|tri| tri.map(|v| 10 * v)).collect();
        assert_eq!(
            t.to_polygon_corners(&values, 4, u32::MAX),
            Some(vec![0, 10, 20, 30])
        );
        assert_eq!(t.to_polygon_corners(&values[..1], 4, 0), None);
        assert_eq!(t.to_polygon_corners(&values, 3, 0), None);
    }
}
