//! The input mesh as the intrinsic triangulation sees it.

use ddg_core::topology::NONE;
use ddg_core::{Geometry, TriMesh};

use crate::error::IntrinsicError;
use crate::formulas::angle_scale;
use crate::point::SurfacePoint;
use crate::surface::Surface;
use crate::trace::{Trace, TraceOptions, trace_query};

/// The input mesh with what tracing over it needs: its connectivity, its edge lengths, the
/// angle sum `Θ` of every vertex and the signpost of every halfedge.
///
/// It is a copy of the caller's `TriMesh` and `Geometry`, made once by
/// [`IntrinsicTriangulation::new`](crate::IntrinsicTriangulation::new) and never changed. The
/// triangulation holds it behind an `Arc`: it has no lifetime parameter (so it can live in a
/// session object or behind a Python handle), and clones of a triangulation share the input.
/// The cost is one more copy of the input's connectivity, which is of the size of the
/// triangulation itself.
///
/// The signposts of the input mesh define the frame of every input vertex: angle zero is the
/// direction of [`TriMesh::outgoing_halfedge`] (at a boundary vertex the boundary halfedge that
/// has the surface on its left), and the signpost of the `n`-th outgoing halfedge
/// counter-clockwise is the sum of the corner angles before it, times
/// [`angle_scale`](crate::formulas::angle_scale).
#[derive(Debug, Clone, PartialEq)]
pub struct InputSurface {
    mesh: TriMesh,
    geometry: Geometry,
    signposts: Vec<f64>,
    angle_sums: Vec<f64>,
}

impl InputSurface {
    /// Copies the mesh and its geometry and computes angle sums and signposts (paper Alg. 4).
    ///
    /// # Errors
    /// [`IntrinsicError::Geometry`] if the geometry was not built for the mesh.
    pub fn new(mesh: &TriMesh, geometry: &Geometry) -> Result<InputSurface, IntrinsicError> {
        geometry.check_mesh(mesh)?;
        let angles = geometry.corner_angles();
        let mut signposts = vec![0.0; mesh.n_halfedges()];
        let mut angle_sums = vec![0.0; mesh.n_vertices()];
        let mut fan: Vec<usize> = Vec::new();
        for (v, angle_sum) in angle_sums.iter_mut().enumerate() {
            let theta: f64 = mesh
                .outgoing_halfedges(v)
                .map(|h| angles[h / 3][h % 3])
                .sum();
            let boundary = mesh.is_boundary_vertex(v);
            let scale = angle_scale(theta, boundary);
            // The first halfedge, whose signpost is 0 (SPEC N9): at a boundary vertex the
            // outgoing halfedge without a twin, where the fan starts; at an interior vertex
            // the outgoing halfedge with the smallest id.
            fan.clear();
            fan.extend(mesh.outgoing_halfedges(v));
            let first = if boundary {
                0
            } else {
                (0..fan.len()).min_by_key(|&i| fan[i]).unwrap_or(0)
            };
            let mut running = 0.0;
            for i in 0..fan.len() {
                let h = fan[(first + i) % fan.len()];
                signposts[h] = running * scale;
                running += angles[h / 3][h % 3];
            }
            *angle_sum = theta;
        }
        Ok(InputSurface {
            mesh: mesh.clone(),
            geometry: geometry.clone(),
            signposts,
            angle_sums,
        })
    }

    /// The connectivity of the input mesh.
    #[inline]
    pub fn mesh(&self) -> &TriMesh {
        &self.mesh
    }

    /// The geometry of the input mesh.
    #[inline]
    pub fn geometry(&self) -> &Geometry {
        &self.geometry
    }

    /// The normalized direction of input halfedge `h` at its tail.
    #[inline]
    pub fn signpost(&self, h: usize) -> f64 {
        self.signposts[h]
    }

    /// The signposts of all input halfedges.
    #[inline]
    pub fn signposts(&self) -> &[f64] {
        &self.signposts
    }

    /// The sum `Θ` of the corner angles at input vertex `v`.
    #[inline]
    pub fn vertex_angle_sum(&self, v: usize) -> f64 {
        self.angle_sums[v]
    }

    /// The angle sums of all input vertices.
    #[inline]
    pub fn vertex_angle_sums(&self) -> &[f64] {
        &self.angle_sums
    }

    /// Follows the straight path that leaves `start` in `direction` over the input mesh for
    /// `length` (paper Sec. 3.2.2, Alg. 2). See [`Trace`] for what is returned and
    /// [`SurfacePoint`] for the frame `direction` refers to.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`], [`IntrinsicError::InvalidDirection`],
    /// [`IntrinsicError::InvalidLength`], [`IntrinsicError::DirectionLeavesFace`] and
    /// [`IntrinsicError::CrossingLimit`].
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

impl Surface for InputSurface {
    const NAME: &'static str = "input mesh";

    #[inline]
    fn n_vertices(&self) -> usize {
        self.mesh.n_vertices()
    }
    #[inline]
    fn n_edges(&self) -> usize {
        self.mesh.n_edges()
    }
    #[inline]
    fn n_faces(&self) -> usize {
        self.mesh.n_faces()
    }
    #[inline]
    fn twin(&self, h: usize) -> u32 {
        self.mesh.twin(h).map_or(NONE, |t| t as u32)
    }
    #[inline]
    fn tail(&self, h: usize) -> u32 {
        self.mesh.from_vertex(h) as u32
    }
    #[inline]
    fn edge(&self, h: usize) -> u32 {
        self.mesh.edge(h) as u32
    }
    #[inline]
    fn length(&self, e: usize) -> f64 {
        self.geometry.edge_lengths()[e]
    }
    #[inline]
    fn outgoing(&self, v: usize) -> u32 {
        self.mesh.outgoing_halfedge(v) as u32
    }
    #[inline]
    fn signpost(&self, h: usize) -> f64 {
        self.signposts[h]
    }
    #[inline]
    fn angle_sum(&self, v: usize) -> f64 {
        self.angle_sums[v]
    }
    #[inline]
    fn is_forward(&self, h: usize) -> bool {
        self.mesh.from_vertex(h) == self.mesh.edge_vertices(self.mesh.edge(h))[0]
    }
    #[inline]
    fn edge_halfedge(&self, e: usize) -> u32 {
        self.mesh.edge_halfedge(e) as u32
    }
}
