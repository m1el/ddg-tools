//! The compacted coarse mesh and the pointwise map (SPEC §7.1) and the colouring of the coarse
//! faces (SPEC §7.4; paper Sec. 7.1 and 8.1).
//!
//! [`CoarseMesh::compact`] copies what the map and the prolongation operators need out of the
//! working structure: faces as ids `0..|V̄|` of the coarse vertices (with the map to the fine
//! ids), the lengths of the sides of every face, the corner angles and angular coordinates of
//! every halfedge, and every tracked point as [`CoarsePoint::Vertex`] or [`CoarsePoint::Face`]
//! in these ids. The coarse mesh needs no renumbering: the Δ-complex keeps its vertex, edge and
//! face ids dense (a removal moves the last element into the freed id), so the compacted ids
//! are the complex's ids, and `faces()` of [`CoarseMesh`] lists the same faces with fine ids.
//!
//! The prolongation operators are in [`crate::prolongation`].

use crate::error::IemError;
use crate::mesh::{CoarseMesh, Location};

/// Where a tracked point is on the compacted coarse mesh (SPEC §7.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoarsePoint {
    /// The point is the coarse vertex with this compacted id.
    Vertex(u32),
    /// The point is in coarse face `face` with barycentric coordinates `bary` w.r.t. the
    /// corners `CoarseOutput::faces[face]` (`bary[k]` belongs to slot `k`), `≥ 0`, summing to
    /// one.
    Face {
        /// The coarse face.
        face: u32,
        /// Its barycentric coordinates.
        bary: [f64; 3],
    },
}

/// The coarse mesh in compacted ids, with the map from the tracked points.
///
/// Halfedge `3f + k` of face `f` goes from slot `k` to slot `k + 1` (SPEC §2.1). All
/// per-halfedge arrays are indexed `3f + k`.
#[derive(Debug, Clone, PartialEq)]
pub struct CoarseOutput {
    /// The number of fine (input) vertices; points `0..n_fine_vertices` are the fine vertices.
    pub n_fine_vertices: usize,
    /// The fine id of every coarse vertex.
    pub vertices: Vec<u32>,
    /// The coarse id of every fine vertex that is a coarse vertex, `None` if it was removed.
    pub coarse_of_fine: Vec<Option<u32>>,
    /// The faces, as coarse vertex ids (a Δ-face may repeat a vertex).
    pub faces: Vec<[u32; 3]>,
    /// The edge of every halfedge (compacted edge ids).
    pub face_edges: Vec<[u32; 3]>,
    /// The two end points of every edge (coarse ids; the first is the tail of the edge's
    /// primary halfedge).
    pub edges: Vec<[u32; 2]>,
    /// The length of every edge.
    pub edge_lengths: Vec<f64>,
    /// Per face: the length of the side from slot `k` to slot `k + 1`.
    pub face_lengths: Vec<[f64; 3]>,
    /// Per face: the corner angle at slot `k`.
    pub corner_angles: Vec<[f64; 3]>,
    /// Per halfedge: the other halfedge of its edge, `None` on the boundary.
    pub twin: Vec<Option<u32>>,
    /// Per halfedge: the normalized angular coordinate at its tail (SPEC S4).
    pub phi: Vec<f64>,
    /// Per coarse vertex: the sum `Θ` of its corner angles (as cached by the operations).
    pub angle_sum: Vec<f64>,
    /// Per coarse vertex: whether it is on the boundary.
    pub boundary: Vec<bool>,
    /// The location of every tracked point (the fine vertices, then the extra points).
    pub points: Vec<CoarsePoint>,
}

impl CoarseMesh {
    /// The compacted coarse mesh with the location of every tracked point (SPEC §7.1).
    ///
    /// # Errors
    /// [`IemError::RemovalPending`] while a removal is pending (the structure is mid-operation).
    pub fn compact(&self) -> Result<CoarseOutput, IemError> {
        if self.pending.is_some() {
            return Err(IemError::RemovalPending);
        }
        let c = &self.complex;
        let nf = c.n_faces();
        let faces = c.faces().to_vec();
        let face_edges = (0..nf)
            .map(|f| [0, 1, 2].map(|k| c.edge(3 * f + k) as u32))
            .collect();
        let edges = (0..c.n_edges())
            .map(|e| c.edge_vertices(e).map(|v| v as u32))
            .collect();
        let face_lengths = (0..nf)
            .map(|f| crate::mesh::face_lengths(c, &self.length, f))
            .collect();
        let corner_angles = (0..nf)
            .map(|f| [0, 1, 2].map(|k| self.corner_angle(3 * f + k)))
            .collect();
        let twin = (0..c.n_halfedges())
            .map(|h| c.twin(h).map(|t| t as u32))
            .collect();
        let points = self
            .points
            .iter()
            .map(|p| match *p {
                Location::Vertex(fine) => CoarsePoint::Vertex(self.complex_of[fine as usize]),
                Location::Face { face, bary } => CoarsePoint::Face { face, bary },
            })
            .collect();
        Ok(CoarseOutput {
            n_fine_vertices: self.n_input as usize,
            vertices: self.fine_of.clone(),
            coarse_of_fine: (0..self.n_input)
                .map(|v| self.complex_vertex(v).map(|c| c as u32))
                .collect(),
            faces,
            face_edges,
            edges,
            edge_lengths: self.length.clone(),
            face_lengths,
            corner_angles,
            twin,
            phi: self.phi.clone(),
            angle_sum: self.angle_sum.clone(),
            boundary: (0..c.n_vertices())
                .map(|v| c.is_boundary_vertex(v))
                .collect(),
            points,
        })
    }
}

/// The greedy colouring of the coarse faces (SPEC §7.4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceColouring {
    /// The colour index of every coarse face; faces that share a vertex have different colours.
    pub colours: Vec<u32>,
    /// The number of colours used (`max colour + 1`, `0` for no faces).
    pub n_colours: u32,
}

impl CoarseOutput {
    /// The number of coarse vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// The number of coarse faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// The faces with fine vertex ids (what [`CoarseMesh::faces`] returns).
    pub fn faces_fine(&self) -> Vec<[u32; 3]> {
        self.faces
            .iter()
            .map(|f| f.map(|v| self.vertices[v as usize]))
            .collect()
    }

    /// The coarse face that contains tracked point `p` (paper Sec. 8.1): its face, or for a
    /// point that is a coarse vertex the incident face of the lowest id (deterministic).
    /// `None` if `p` is out of range or the vertex has no face.
    pub fn face_of_point(&self, p: usize) -> Option<u32> {
        match *self.points.get(p)? {
            CoarsePoint::Face { face, .. } => Some(face),
            CoarsePoint::Vertex(v) => self
                .faces
                .iter()
                .position(|f| f.contains(&v))
                .map(|f| f as u32),
        }
    }

    /// [`face_of_point`](Self::face_of_point) for every tracked point, with one pass over the
    /// faces for the points that are vertices.
    pub fn face_of_points(&self) -> Vec<Option<u32>> {
        let mut first = vec![None; self.vertices.len()];
        for (f, corners) in self.faces.iter().enumerate() {
            for &v in corners {
                first[v as usize].get_or_insert(f as u32);
            }
        }
        self.points
            .iter()
            .map(|p| match *p {
                CoarsePoint::Face { face, .. } => Some(face),
                CoarsePoint::Vertex(v) => first[v as usize],
            })
            .collect()
    }

    /// The faces incident to every coarse vertex (CSR: `ptr[v]..ptr[v + 1]` indexes `list`;
    /// a face that uses a vertex twice is listed twice). Internal helper of the colouring and
    /// the prolongation.
    pub(crate) fn vertex_slots(&self) -> (Vec<usize>, Vec<(u32, u8)>) {
        let nv = self.vertices.len();
        let mut ptr = vec![0usize; nv + 1];
        for f in &self.faces {
            for &v in f {
                ptr[v as usize + 1] += 1;
            }
        }
        for v in 0..nv {
            ptr[v + 1] += ptr[v];
        }
        let mut next = ptr.clone();
        let mut list = vec![(0u32, 0u8); ptr[nv]];
        for (f, corners) in self.faces.iter().enumerate() {
            for (k, &v) in corners.iter().enumerate() {
                list[next[v as usize]] = (f as u32, k as u8);
                next[v as usize] += 1;
            }
        }
        (ptr, list)
    }

    /// A deterministic greedy colouring of the coarse faces such that faces that share a vertex
    /// (the faces "adjacent through a vertex", Sharp et al. 2019, paper Sec. 8.1) get different
    /// colours: the faces in id order, each taking the smallest colour not used by an
    /// already-coloured face around its corners. There is no palette limit; the caller maps
    /// the indices to colours (the reference's 24-colour overflow, R40, cannot occur).
    pub fn colour_faces(&self) -> FaceColouring {
        let (ptr, list) = self.vertex_slots();
        let nf = self.faces.len();
        const UNSET: u32 = u32::MAX;
        let mut colours = vec![UNSET; nf];
        let mut used: Vec<u32> = Vec::new();
        let mut n_colours = 0u32;
        for f in 0..nf {
            used.clear();
            for &v in &self.faces[f] {
                for &(g, _) in &list[ptr[v as usize]..ptr[v as usize + 1]] {
                    let c = colours[g as usize];
                    if c != UNSET {
                        used.push(c);
                    }
                }
            }
            used.sort_unstable();
            used.dedup();
            let mut colour = 0u32;
            for &u in &used {
                if u == colour {
                    colour += 1;
                } else if u > colour {
                    break;
                }
            }
            colours[f] = colour;
            n_colours = n_colours.max(colour + 1);
        }
        FaceColouring { colours, n_colours }
    }

    /// The colour of the coarse face of every tracked point (`None` for a point without a
    /// face), by [`face_of_points`](Self::face_of_points) and [`colour_faces`](Self::colour_faces):
    /// what the Blender view paints the fine faces with (paper Sec. 8.1).
    pub fn point_colours(&self) -> Vec<Option<u32>> {
        let colouring = self.colour_faces();
        self.face_of_points()
            .into_iter()
            .map(|f| f.map(|f| colouring.colours[f as usize]))
            .collect()
    }
}
