//! Triangle-soup obstacles with a BVH and point–triangle distance (`SPEC.md` §5.3, D19).
//!
//! An [`Obstacle`] is a rigid triangle soup: no manifold, orientation or closedness
//! requirement (polygons of an OBJ file are fan-triangulated, nothing else is checked except
//! finite coordinates, valid indices and non-degenerate triangles). It serves
//!
//! * the surface potential of [`crate::penalties::surface`]: [`Obstacle::potential_at`]
//!   evaluates `φ(x) = Σ_t A_t δ(x, t)^{−p}` with the exact point–triangle distance `δ` in the
//!   near field and centroid lumping of far BVH nodes (the BVH is the generic
//!   [`crate::spatial::bvh::Bvh`] over the triangle centroids with the areas as masses);
//! * the continuous collision detection, through [`CollisionObstacle`].
//!
//! Closest points follow Ericson, *Real-Time Collision Detection* (2004), §5.1.5
//! [formula: rc-closest-points]. The distance to a triangle (a closed convex set) is `C¹`
//! outside it with `∇δ = (x − π_t(x))/δ` [formula: rc-distance-gradient].
//!
//! The authors' obstacle code (`src/obstacles/mesh_obstacle.cpp`, MIT, Copyright (c) 2019
//! Christopher Yu) samples the mesh vertices with a BVH of point masses; this one puts the
//! barrier on the faces (D19).
//!
//! **Determinism.** A query depends on its point only (fixed traversal order); the callers
//! reduce per-vertex results in index order, so results do not depend on the thread count.

use std::path::Path;

use crate::ccd::CollisionObstacle;
use crate::error::RcError;
use crate::kernel::v3::{add, dot, norm, scale, sub};
use crate::spatial::aabb::Aabb;
use crate::spatial::bvh::{Bvh, Element};

/// The default BVH leaf size of an obstacle.
pub const OBSTACLE_LEAF_SIZE: usize = 4;

/// The closest point of the triangle `a b c` to `p` (Ericson 2004, §5.1.5). In exact arithmetic,
/// for a non-degenerate triangle, it is the closest point (`DDG.RC.closestPtTriRC_eq`).
///
/// [formula: rc-closest-points]
pub fn closest_point_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return add(a, scale(ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return add(a, scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    add(add(a, scale(ab, vb * denom)), scale(ac, vc * denom))
}

/// Squared distance from `p` to the box.
fn box_distance_sq(b: &Aabb<3>, p: [f64; 3]) -> f64 {
    let mut s = 0.0;
    for k in 0..3 {
        let d = (b.lo[k] - p[k]).max(0.0).max(p[k] - b.hi[k]);
        s += d * d;
    }
    s
}

#[derive(Debug, Clone)]
struct TriElem {
    corners: [[f64; 3]; 3],
    centroid: [f64; 3],
    area: f64,
}

impl Element<3> for TriElem {
    fn point(&self) -> [f64; 3] {
        self.centroid
    }
    fn bounds(&self) -> Aabb<3> {
        let mut b = Aabb::empty();
        for c in self.corners {
            b.grow(c);
        }
        b
    }
    fn mass(&self) -> f64 {
        self.area
    }
}

/// The result of a potential query at one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PotentialValue {
    /// `φ(x) = Σ_t A_t δ^{−p}` (`+∞` if `x` lies on a triangle).
    pub value: f64,
    /// `∇φ(x)` (zero if not requested; NaN if `x` lies on a triangle).
    pub gradient: [f64; 3],
    /// Number of exact point–triangle evaluations (near field).
    pub near_triangles: u32,
    /// Number of lumped BVH nodes (far field).
    pub far_nodes: u32,
}

/// A rigid triangle soup with a BVH (`SPEC.md` §5.3, D19).
#[derive(Debug, Clone)]
pub struct Obstacle {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    face_of_triangle: Vec<u32>,
    corners: Vec<[[f64; 3]; 3]>,
    areas: Vec<f64>,
    bvh: Bvh<3>,
}

impl CollisionObstacle for Obstacle {
    fn triangles(&self) -> &[[[f64; 3]; 3]] {
        &self.corners
    }
}

fn bad(name: &'static str, value: String, requirement: &'static str) -> RcError {
    RcError::InvalidOption {
        name,
        value,
        requirement,
    }
}

impl Obstacle {
    /// An obstacle from indexed triangles.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for an empty triangle list;
    /// [`RcError::NonFinitePosition`] for a NaN or infinite coordinate;
    /// [`RcError::ObstacleTriangleOutOfRange`] for an index `>= vertices.len()`;
    /// [`RcError::DegenerateObstacleTriangle`] for a triangle of (numerically) zero area, i.e.
    /// area `<= 1e-12 ×` the square of its longest edge.
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> Result<Obstacle, RcError> {
        let faces: Vec<u32> = (0..triangles.len() as u32).collect();
        Obstacle::with_faces(vertices, triangles, faces)
    }

    fn with_faces(
        vertices: Vec<[f64; 3]>,
        triangles: Vec<[u32; 3]>,
        face_of_triangle: Vec<u32>,
    ) -> Result<Obstacle, RcError> {
        if triangles.is_empty() {
            return Err(bad(
                "obstacle.triangles",
                "[]".into(),
                "at least one triangle",
            ));
        }
        if triangles.len() > u32::MAX as usize {
            return Err(bad(
                "obstacle.triangles",
                triangles.len().to_string(),
                "at most 2^32 - 1 triangles",
            ));
        }
        for (v, p) in vertices.iter().enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                return Err(RcError::NonFinitePosition { vertex: v as u32 });
            }
        }
        let mut corners = Vec::with_capacity(triangles.len());
        let mut areas = Vec::with_capacity(triangles.len());
        let mut elems = Vec::with_capacity(triangles.len());
        for (t, tri) in triangles.iter().enumerate() {
            for &v in tri {
                if v as usize >= vertices.len() {
                    return Err(RcError::ObstacleTriangleOutOfRange {
                        triangle: t,
                        vertex: v,
                        n_vertices: vertices.len(),
                    });
                }
            }
            let c = tri.map(|v| vertices[v as usize]);
            let (e1, e2, e3) = (sub(c[1], c[0]), sub(c[2], c[0]), sub(c[2], c[1]));
            let n = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let area = 0.5 * norm(n);
            let longest = dot(e1, e1).max(dot(e2, e2)).max(dot(e3, e3));
            if area.is_nan() || area <= 1e-12 * longest {
                return Err(RcError::DegenerateObstacleTriangle { triangle: t });
            }
            let centroid = scale(add(add(c[0], c[1]), c[2]), 1.0 / 3.0);
            corners.push(c);
            areas.push(area);
            elems.push(TriElem {
                corners: c,
                centroid,
                area,
            });
        }
        let bvh = Bvh::build(&elems, OBSTACLE_LEAF_SIZE)?;
        Ok(Obstacle {
            vertices,
            triangles,
            face_of_triangle,
            corners,
            areas,
            bvh,
        })
    }

    /// An obstacle from a triangle soup (corner coordinates).
    ///
    /// # Errors
    /// As [`Obstacle::new`].
    pub fn from_soup(soup: &[[[f64; 3]; 3]]) -> Result<Obstacle, RcError> {
        let vertices: Vec<[f64; 3]> = soup.iter().flatten().copied().collect();
        let triangles = (0..soup.len() as u32)
            .map(|t| [3 * t, 3 * t + 1, 3 * t + 2])
            .collect();
        Obstacle::new(vertices, triangles)
    }

    /// An obstacle from an OBJ mesh of `ddg_core`: every polygon `(c₀, …, c_{k−1})` is
    /// fan-triangulated as `(c₀, c_i, c_{i+1})`; no manifold, orientation or closedness check.
    ///
    /// # Errors
    /// As [`Obstacle::new`] (`DegenerateObstacleTriangle` carries the triangle's index in the
    /// fan-triangulated list, see [`Obstacle::face_of_triangle`] for the polygon).
    pub fn from_obj_mesh(mesh: &ddg_core::io::ObjMesh) -> Result<Obstacle, RcError> {
        let mut tris = Vec::new();
        let mut faces = Vec::new();
        for f in 0..mesh.n_faces() {
            let r = mesh.face_corners(f);
            let c = &mesh.corner_verts[r];
            for i in 1..c.len().saturating_sub(1) {
                tris.push([c[0], c[i], c[i + 1]]);
                faces.push(f as u32);
            }
        }
        Obstacle::with_faces(mesh.positions.clone(), tris, faces)
    }

    /// Parses OBJ text (`v` and `f` statements; see `ddg_core::io::parse_obj`).
    ///
    /// # Errors
    /// [`RcError::Parse`] for malformed text (with the line number); then as
    /// [`Obstacle::from_obj_mesh`]. A file without faces is an [`RcError::InvalidOption`].
    pub fn parse_obj(text: &str) -> Result<Obstacle, RcError> {
        let mesh = ddg_core::io::parse_obj(text).map_err(convert_file_error)?;
        Obstacle::from_obj_mesh(&mesh)
    }

    /// Reads an OBJ file with `ddg_core::io::read_obj`.
    ///
    /// # Errors
    /// [`RcError::Io`], [`RcError::Parse`] (the message names the file), then as
    /// [`Obstacle::from_obj_mesh`].
    pub fn read_obj(path: impl AsRef<Path>) -> Result<Obstacle, RcError> {
        let mesh = ddg_core::io::read_obj(path).map_err(convert_file_error)?;
        Obstacle::from_obj_mesh(&mesh)
    }

    /// The vertices.
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices
    }

    /// The triangles (vertex indices; fan-triangulated polygons).
    pub fn triangle_indices(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// For every triangle the polygon of the OBJ file it came from (the triangle index itself
    /// for [`Obstacle::new`]).
    pub fn face_of_triangle(&self) -> &[u32] {
        &self.face_of_triangle
    }

    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// The triangle areas.
    pub fn areas(&self) -> &[f64] {
        &self.areas
    }

    /// The BVH over the triangle centroids (masses: the areas).
    pub fn bvh(&self) -> &Bvh<3> {
        &self.bvh
    }

    /// The distance `δ(x, t)` and the closest point of triangle `t`.
    pub fn triangle_distance(&self, x: [f64; 3], t: usize) -> (f64, [f64; 3]) {
        let c = &self.corners[t];
        let q = closest_point_on_triangle(x, c[0], c[1], c[2]);
        (norm(sub(x, q)), q)
    }

    /// The distance from `x` to the soup, by branch and bound over the BVH (ties: lower
    /// triangle index). Returns `(distance, triangle)`.
    pub fn distance_to(&self, x: [f64; 3]) -> (f64, u32) {
        let nodes = self.bvh.nodes();
        let mut best = (f64::INFINITY, u32::MAX);
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let nd = &nodes[i];
            if box_distance_sq(&nd.bbox, x).sqrt() > best.0 {
                continue;
            }
            if nd.is_leaf() {
                for &t in self.bvh.elements_of(i) {
                    let (d, _) = self.triangle_distance(x, t as usize);
                    if d < best.0 || (d == best.0 && t < best.1) {
                        best = (d, t);
                    }
                }
            } else {
                let (l, r) = (i + 1, nd.right as usize);
                let dl = box_distance_sq(&nodes[l].bbox, x);
                let dr = box_distance_sq(&nodes[r].bbox, x);
                // the nearer child is popped first
                if dl <= dr {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            }
        }
        best
    }

    /// The distance from `x` to the soup by brute force over all triangles.
    pub fn distance_to_brute(&self, x: [f64; 3]) -> f64 {
        (0..self.n_triangles())
            .map(|t| self.triangle_distance(x, t).0)
            .fold(f64::INFINITY, f64::min)
    }

    /// `φ(x) = Σ_t A_t δ(x, t)^{−p}` and optionally `∇φ`, by the BVH: a node is used as one
    /// lumped term `A_N |x − c_N|^{−p}` (`c_N` the area-weighted centroid) when `x` lies
    /// outside its box and `r_N <= θ |x − c_N|` (`r_N` the half diagonal of the box); otherwise
    /// its triangles are visited exactly. `θ = 0` is exact. At a point on a triangle the value
    /// is `+∞` and the gradient NaN.
    ///
    /// [formula: rc-penalty-surface] [formula: rc-distance-gradient]
    pub fn potential_at(
        &self,
        x: [f64; 3],
        p: f64,
        theta: f64,
        want_gradient: bool,
    ) -> PotentialValue {
        let nodes = self.bvh.nodes();
        let mut out = PotentialValue {
            value: 0.0,
            gradient: [0.0; 3],
            near_triangles: 0,
            far_nodes: 0,
        };
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let nd = &nodes[i];
            if theta > 0.0 && !nd.bbox.contains_point(x) {
                let c = nd.centroid();
                let d = sub(x, c);
                let r = norm(d);
                if nd.bbox.half_diagonal() <= theta * r {
                    out.value += nd.mass * r.powf(-p);
                    if want_gradient {
                        let f = -p * nd.mass * r.powf(-p - 2.0);
                        out.gradient = add(out.gradient, scale(d, f));
                    }
                    out.far_nodes += 1;
                    continue;
                }
            }
            if nd.is_leaf() {
                for &t in self.bvh.elements_of(i) {
                    self.add_triangle(&mut out, x, t as usize, p, want_gradient);
                }
            } else {
                stack.push(nd.right as usize);
                stack.push(i + 1);
            }
        }
        out
    }

    /// As [`Obstacle::potential_at`] with every triangle exact (reference for tests).
    pub fn potential_at_brute(&self, x: [f64; 3], p: f64, want_gradient: bool) -> PotentialValue {
        let mut out = PotentialValue {
            value: 0.0,
            gradient: [0.0; 3],
            near_triangles: 0,
            far_nodes: 0,
        };
        for t in 0..self.n_triangles() {
            self.add_triangle(&mut out, x, t, p, want_gradient);
        }
        out
    }

    fn add_triangle(&self, out: &mut PotentialValue, x: [f64; 3], t: usize, p: f64, grad: bool) {
        let (d, q) = self.triangle_distance(x, t);
        let a = self.areas[t];
        out.near_triangles += 1;
        if d == 0.0 {
            out.value = f64::INFINITY;
            out.gradient = [f64::NAN; 3];
            return;
        }
        out.value += a * d.powf(-p);
        if grad {
            // ∇ δ^{-p} = -p δ^{-p-1} (x - q)/δ
            let f = -p * a * d.powf(-p - 2.0);
            out.gradient = add(out.gradient, scale(sub(x, q), f));
        }
    }
}

fn convert_file_error(e: ddg_core::io::FileError) -> RcError {
    use ddg_core::io::FileError;
    match e {
        FileError::Io { path, message } => RcError::Io { path, message },
        FileError::Parse {
            path,
            line,
            message,
        } => RcError::Parse {
            line,
            message: match path {
                Some(p) => format!("{p}: {message}"),
                None => message,
            },
        },
        FileError::InvalidInput { message } => RcError::Parse { line: 0, message },
    }
}
