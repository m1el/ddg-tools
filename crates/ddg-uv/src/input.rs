//! The input of the pipeline and the checks on its arrays.

use std::f64::consts::TAU;

use ddg_core::io::Cone;

use crate::error::UvError;

/// A polygon mesh, its triangulation and the user's steering. Everything is borrowed; nothing
/// is modified.
///
/// * Polygon `p` owns the corners `polygon_offsets[p] .. polygon_offsets[p + 1]`, in winding
///   order; `corner_vertices[c]` is the vertex of corner `c`.
/// * `triangles` triangulates every polygon: each triangle is three **corner** indices of one
///   polygon, wound like the polygon; every side of the polygon is in exactly one triangle and
///   every diagonal in exactly two. A polygon with `n` corners has `n − 2` triangles. The
///   triangles may come in any order. For a triangle mesh this is `[3p, 3p + 1, 3p + 2]`.
///   Blender's `mesh.loop_triangles` satisfies it. The pipeline never chooses a triangulation
///   itself.
/// * `seam_edges`, `pinned_vertices` and `forbidden_vertices` are sets: order, orientation and
///   repetitions do not matter.
///
/// [`MeshInput::validate`] checks all of that; [`crate::unwrap`] calls it.
#[derive(Debug, Clone, Copy)]
pub struct MeshInput<'a> {
    /// Vertex positions.
    pub positions: &'a [[f64; 3]],
    /// The vertex of every polygon corner.
    pub corner_vertices: &'a [u32],
    /// Corner offsets of the polygons: `n_polygons + 1` entries, the first 0, the last the
    /// number of corners.
    pub polygon_offsets: &'a [u32],
    /// The triangulation, as triples of polygon-corner indices.
    pub triangles: &'a [[u32; 3]],
    /// Seams as vertex pairs. The surface is cut along them before anything else.
    pub seam_edges: &'a [[u32; 2]],
    /// Vertices where a cone costs nothing; its angle is free (a free cone site, paper
    /// Sec. 8.1).
    pub pinned_vertices: &'a [u32],
    /// Cones whose vertex and angle are fixed.
    pub prescribed_cones: &'a [Cone],
    /// Vertices where no cone may be placed.
    pub forbidden_vertices: &'a [u32],
}

impl<'a> MeshInput<'a> {
    /// A mesh without seams and steering.
    pub fn new(
        positions: &'a [[f64; 3]],
        corner_vertices: &'a [u32],
        polygon_offsets: &'a [u32],
        triangles: &'a [[u32; 3]],
    ) -> MeshInput<'a> {
        MeshInput {
            positions,
            corner_vertices,
            polygon_offsets,
            triangles,
            seam_edges: &[],
            pinned_vertices: &[],
            prescribed_cones: &[],
            forbidden_vertices: &[],
        }
    }

    /// Number of polygons (0 if `polygon_offsets` is empty).
    pub fn n_polygons(&self) -> usize {
        self.polygon_offsets.len().saturating_sub(1)
    }

    /// Number of polygon corners.
    pub fn n_corners(&self) -> usize {
        self.corner_vertices.len()
    }

    /// Checks the arrays: lengths, index ranges, finite numbers, that the triangles
    /// triangulate the polygons, and that the steering is not contradictory.
    ///
    /// Not checked here, because it needs the connectivity: whether the mesh is manifold and
    /// orientable, whether a polygon repeats a vertex, whether a seam is an edge of the mesh
    /// (see the stages of [`crate::unwrap`]).
    ///
    /// # Errors
    /// The input variants of [`UvError`] (`kind() == ErrorKind::Input`).
    pub fn validate(&self) -> Result<(), UvError> {
        const LIMIT: usize = u32::MAX as usize;
        let n_vertices = self.positions.len();
        let n_corners = self.corner_vertices.len();
        if n_vertices >= LIMIT {
            return Err(UvError::TooLarge {
                what: "positions",
                count: n_vertices,
            });
        }
        if n_corners >= LIMIT {
            return Err(UvError::TooLarge {
                what: "corner_vertices",
                count: n_corners,
            });
        }
        for (v, p) in self.positions.iter().enumerate() {
            if let Some(&x) = p.iter().find(|x| !x.is_finite()) {
                return Err(UvError::NonFinite {
                    what: "positions",
                    position: v,
                    value: x,
                });
            }
        }

        // Polygons.
        let offsets = self.polygon_offsets;
        if offsets.len() < 2 {
            return Err(UvError::NoPolygons);
        }
        if offsets[0] != 0 {
            return Err(UvError::PolygonOffsets {
                position: 0,
                value: offsets[0],
                requirement: "the first face must start at corner 0".to_string(),
            });
        }
        let n_polygons = offsets.len() - 1;
        for p in 0..n_polygons {
            let (start, end) = (offsets[p], offsets[p + 1]);
            if end < start {
                return Err(UvError::PolygonOffsets {
                    position: p + 1,
                    value: end,
                    requirement: format!(
                        "offsets must not decrease, but the previous one is {start}"
                    ),
                });
            }
            if end - start < 3 {
                return Err(UvError::PolygonTooSmall {
                    polygon: p as u32,
                    n_corners: (end - start) as usize,
                });
            }
        }
        if offsets[n_polygons] as usize != n_corners {
            return Err(UvError::PolygonOffsets {
                position: n_polygons,
                value: offsets[n_polygons],
                requirement: format!(
                    "the last offset must be the number of corners, {n_corners}"
                ),
            });
        }
        check_range(
            "corner_vertices",
            self.corner_vertices.iter().copied(),
            1,
            n_vertices,
            "vertices",
        )?;

        // Triangles.
        check_range(
            "triangles",
            self.triangles.iter().flatten().copied(),
            3,
            n_corners,
            "polygon corners",
        )?;
        let expected = n_corners - 2 * n_polygons;
        if self.triangles.len() != expected {
            return Err(UvError::TriangleCount {
                got: self.triangles.len(),
                expected,
            });
        }
        check_triangulation(offsets, self.triangles)?;

        // Seams.
        check_range(
            "seam_edges",
            self.seam_edges.iter().flatten().copied(),
            2,
            n_vertices,
            "vertices",
        )?;
        if let Some(position) = self.seam_edges.iter().position(|e| e[0] == e[1]) {
            return Err(UvError::DegenerateSeamEdge {
                position,
                vertex: self.seam_edges[position][0],
            });
        }

        // Steering.
        check_range(
            "pinned_vertices",
            self.pinned_vertices.iter().copied(),
            1,
            n_vertices,
            "vertices",
        )?;
        check_range(
            "forbidden_vertices",
            self.forbidden_vertices.iter().copied(),
            1,
            n_vertices,
            "vertices",
        )?;
        check_range(
            "prescribed_cones",
            self.prescribed_cones.iter().map(|c| c.vertex),
            1,
            n_vertices,
            "vertices",
        )?;
        for (position, cone) in self.prescribed_cones.iter().enumerate() {
            if !cone.angle.is_finite() {
                return Err(UvError::NonFinite {
                    what: "prescribed_cones",
                    position,
                    value: cone.angle,
                });
            }
            if cone.angle >= TAU {
                return Err(UvError::PrescribedAngleTooLarge {
                    vertex: cone.vertex,
                    angle: cone.angle,
                });
            }
        }
        // 1 = pinned, 2 = prescribed, 4 = forbidden.
        let mut role = vec![0u8; n_vertices];
        for cone in self.prescribed_cones {
            if role[cone.vertex as usize] & 2 != 0 {
                return Err(UvError::DuplicatePrescribed {
                    vertex: cone.vertex,
                });
            }
            role[cone.vertex as usize] |= 2;
        }
        for &v in self.pinned_vertices {
            role[v as usize] |= 1;
        }
        for &v in self.forbidden_vertices {
            role[v as usize] |= 4;
        }
        if let Some(v) = role.iter().position(|r| r.count_ones() > 1) {
            let names: Vec<&'static str> = [(1, "pinned"), (2, "prescribed"), (4, "forbidden")]
                .iter()
                .filter(|(bit, _)| role[v] & bit != 0)
                .map(|&(_, name)| name)
                .collect();
            return Err(UvError::ConflictingSteering {
                vertex: v as u32,
                first: names[0],
                second: names[1],
            });
        }
        Ok(())
    }
}

/// Checks that every value is below `count`. `width` is the number of values per row; the
/// reported position is the row.
fn check_range(
    what: &'static str,
    values: impl Iterator<Item = u32>,
    width: usize,
    count: usize,
    of: &'static str,
) -> Result<(), UvError> {
    for (i, value) in values.enumerate() {
        if value as usize >= count {
            return Err(UvError::IndexOutOfRange {
                what,
                position: i / width,
                value,
                count,
                of,
            });
        }
    }
    Ok(())
}

/// Checks that the triangles of every polygon triangulate it: `n − 2` triangles of corners of
/// that polygon, wound like it, every side of the polygon in exactly one triangle and every
/// diagonal in exactly two (once in each direction). Offsets and corner indices are already
/// range-checked.
fn check_triangulation(offsets: &[u32], triangles: &[[u32; 3]]) -> Result<(), UvError> {
    let n_polygons = offsets.len() - 1;
    let n_corners = offsets[n_polygons] as usize;
    let mut corner_polygon = vec![0u32; n_corners];
    for p in 0..n_polygons {
        corner_polygon[offsets[p] as usize..offsets[p + 1] as usize].fill(p as u32);
    }
    let bad = |polygon: u32, problem: String| UvError::BadTriangulation { polygon, problem };
    let next = |c: u32, p: u32| {
        if c + 1 == offsets[p as usize + 1] {
            offsets[p as usize]
        } else {
            c + 1
        }
    };
    let mut triangle_count = vec![0u32; n_polygons];
    // Per corner: how many triangles contain the polygon side that starts there.
    let mut side_count = vec![0u8; n_corners];
    // Diagonals as (lower corner, higher corner, runs from lower to higher).
    let mut diagonals: Vec<(u32, u32, bool)> = Vec::new();
    for (t, &[a, b, c]) in triangles.iter().enumerate() {
        let p = corner_polygon[a as usize];
        let (pb, pc) = (corner_polygon[b as usize], corner_polygon[c as usize]);
        if pb != p || pc != p {
            return Err(bad(
                p,
                format!(
                    "triangle {t} = corners ({a}, {b}, {c}) spans the faces {p}, {pb} and {pc}"
                ),
            ));
        }
        if a == b || b == c || a == c {
            return Err(bad(
                p,
                format!("triangle {t} = corners ({a}, {b}, {c}) repeats a corner"),
            ));
        }
        triangle_count[p as usize] += 1;
        for (from, to) in [(a, b), (b, c), (c, a)] {
            if next(from, p) == to {
                side_count[from as usize] = side_count[from as usize].saturating_add(1);
            } else if next(to, p) == from {
                return Err(bad(
                    p,
                    format!(
                        "triangle {t} = corners ({a}, {b}, {c}) is wound against the face"
                    ),
                ));
            } else {
                diagonals.push((from.min(to), from.max(to), from < to));
            }
        }
    }
    for (p, &count) in triangle_count.iter().enumerate() {
        let size = (offsets[p + 1] - offsets[p]) as usize;
        if count as usize != size - 2 {
            return Err(bad(
                p as u32,
                format!(
                    "it has {size} corners and needs {} triangles, but {count} were given",
                    size - 2
                ),
            ));
        }
    }
    if let Some(corner) = side_count.iter().position(|&n| n != 1) {
        return Err(bad(
            corner_polygon[corner],
            format!(
                "the side that starts at corner {corner} is in {} triangles instead of 1",
                side_count[corner]
            ),
        ));
    }
    // After sorting, the two halves of a diagonal are adjacent, the one running from the
    // higher to the lower corner first.
    diagonals.sort_unstable();
    let mut i = 0;
    while i < diagonals.len() {
        let (lo, hi, forward) = diagonals[i];
        let paired = !forward
            && i + 1 < diagonals.len()
            && diagonals[i + 1] == (lo, hi, true)
            && (i + 2 == diagonals.len() || (diagonals[i + 2].0, diagonals[i + 2].1) != (lo, hi));
        if !paired {
            return Err(bad(
                corner_polygon[lo as usize],
                format!(
                    "the diagonal between corners {lo} and {hi} is not shared by exactly two \
                     triangles, one on each side"
                ),
            ));
        }
        i += 2;
    }
    Ok(())
}

/// An owned polygon mesh with a triangulation: what the generators of [`crate::primitives`]
/// return, and a convenient way to build a [`MeshInput`].
#[derive(Debug, Clone, PartialEq)]
pub struct PolygonMesh {
    /// Vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// The vertex of every polygon corner.
    pub corner_vertices: Vec<u32>,
    /// Corner offsets of the polygons, `n_polygons + 1` entries.
    pub polygon_offsets: Vec<u32>,
    /// The triangulation as triples of polygon-corner indices.
    pub triangles: Vec<[u32; 3]>,
}

impl PolygonMesh {
    /// A triangle mesh: every triangle is a polygon and its own triangulation.
    pub fn from_triangles(positions: Vec<[f64; 3]>, faces: &[[u32; 3]]) -> PolygonMesh {
        let n = faces.len() as u32;
        PolygonMesh {
            positions,
            corner_vertices: faces.iter().flatten().copied().collect(),
            polygon_offsets: (0..=n).map(|f| 3 * f).collect(),
            triangles: (0..n).map(|f| [3 * f, 3 * f + 1, 3 * f + 2]).collect(),
        }
    }

    /// A polygon mesh with every polygon split into a fan around its first corner. A fan is a
    /// valid triangulation of any polygon in the combinatorial sense that [`MeshInput`]
    /// requires; its triangles have positive area if the polygon is convex (more generally,
    /// star-shaped as seen from its first corner). This is not Blender's tessellation.
    pub fn from_polygons(positions: Vec<[f64; 3]>, polygons: &[Vec<u32>]) -> PolygonMesh {
        let mut corner_vertices = Vec::new();
        let mut polygon_offsets = vec![0u32];
        for polygon in polygons {
            corner_vertices.extend_from_slice(polygon);
            polygon_offsets.push(corner_vertices.len() as u32);
        }
        let triangles = fan_triangulation(&polygon_offsets);
        PolygonMesh {
            positions,
            corner_vertices,
            polygon_offsets,
            triangles,
        }
    }

    /// The mesh as pipeline input, without seams and steering.
    pub fn input(&self) -> MeshInput<'_> {
        MeshInput::new(
            &self.positions,
            &self.corner_vertices,
            &self.polygon_offsets,
            &self.triangles,
        )
    }

    /// Number of polygons.
    pub fn n_polygons(&self) -> usize {
        self.polygon_offsets.len().saturating_sub(1)
    }

    /// The corners of polygon `p`.
    pub fn polygon_corners(&self, p: usize) -> std::ops::Range<usize> {
        self.polygon_offsets[p] as usize..self.polygon_offsets[p + 1] as usize
    }

    /// The vertices of polygon `p`, in winding order.
    pub fn polygon(&self, p: usize) -> &[u32] {
        &self.corner_vertices[self.polygon_corners(p)]
    }

    /// The polygon edges as sorted, unique vertex pairs `[min, max]`.
    pub fn edges(&self) -> Vec<[u32; 2]> {
        let mut edges = Vec::with_capacity(self.corner_vertices.len());
        for p in 0..self.n_polygons() {
            let polygon = self.polygon(p);
            for (i, &a) in polygon.iter().enumerate() {
                let b = polygon[(i + 1) % polygon.len()];
                edges.push([a.min(b), a.max(b)]);
            }
        }
        edges.sort_unstable();
        edges.dedup();
        edges
    }
}

/// The fan triangulation of polygons given by their corner offsets: polygon corners
/// `s, s + 1, …, e − 1` give the triangles `(s, i, i + 1)`.
pub fn fan_triangulation(polygon_offsets: &[u32]) -> Vec<[u32; 3]> {
    let mut triangles = Vec::new();
    for w in polygon_offsets.windows(2) {
        let (start, end) = (w[0], w[1]);
        for i in start + 1..end.saturating_sub(1) {
            triangles.push([start, i, i + 1]);
        }
    }
    triangles
}
