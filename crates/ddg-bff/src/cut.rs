//! Cutting a mesh along edges.
//!
//! [`cut_mesh`] splits a [`TriMesh`] along a set of interior edges. Every vertex is replaced by
//! one copy per *sector*: a maximal run of consecutive faces around the vertex that is not
//! interrupted by a cut edge (or by the boundary). A cut edge becomes two boundary edges, one
//! on each side.
//!
//! # Index maps
//!
//! * **Faces and corners keep their ids.** Face `f` of the cut mesh is face `f` of the original
//!   mesh with the same corner slots, so corner (and halfedge) `3 f + k` means the same thing in
//!   both meshes. Per-corner results on the cut mesh are per-corner results on the original.
//! * **Vertices**: [`CutMesh::vertex_to_original`] maps each vertex of the cut mesh to the
//!   vertex it is a copy of. The first sector of original vertex `v` keeps the id `v`;
//!   additional sectors get the ids `n_original_vertices, n_original_vertices + 1, …` in order
//!   of ascending `v` and, for one `v`, counter-clockwise.
//! * **Edges**: [`CutMesh::edge_to_original`]. The two sides of a cut edge map to the same
//!   original edge.

use ddg_core::topology::{NONE, TriMesh};

use crate::error::BffError;

/// A mesh cut along edges, with the maps back to the mesh it was cut from.
///
/// See the [module documentation](self) for the index conventions.
#[derive(Debug, Clone, PartialEq)]
pub struct CutMesh {
    /// The cut mesh. It can have more connected components than the original.
    pub mesh: TriMesh,
    /// For each vertex of `mesh`, the original vertex it is a copy of.
    /// `vertex_to_original[v] == v` for `v < n_original_vertices`.
    pub vertex_to_original: Vec<u32>,
    /// For each edge of `mesh`, the original edge.
    pub edge_to_original: Vec<u32>,
    /// The original interior edges that were cut, ascending and without repetitions.
    pub cut_edges: Vec<u32>,
    /// Number of distinct requested edges that are boundary edges of the original mesh.
    /// Cutting a boundary edge changes nothing; such requests are counted here and ignored.
    pub ignored_boundary_edges: usize,
    /// Number of vertices of the original mesh.
    pub n_original_vertices: usize,
}

impl CutMesh {
    /// The original face of face `f` of the cut mesh. Faces keep their ids, so this is `f`.
    #[inline]
    pub fn face_to_original(&self, f: usize) -> usize {
        f
    }

    /// The original corner of corner `c = 3 f + k` of the cut mesh. Corners keep their ids, so
    /// this is `c`.
    #[inline]
    pub fn corner_to_original(&self, c: usize) -> usize {
        c
    }

    /// Number of vertex copies that cutting added: `n_vertices(cut) − n_vertices(original)`.
    #[inline]
    pub fn n_added_vertices(&self) -> usize {
        self.mesh.n_vertices() - self.n_original_vertices
    }

    /// For each original vertex, how many copies it has in the cut mesh (at least 1).
    pub fn copies_per_vertex(&self) -> Vec<u32> {
        let mut copies = vec![0u32; self.n_original_vertices];
        for &v in &self.vertex_to_original {
            copies[v as usize] += 1;
        }
        copies
    }

    /// Copies a per-vertex array of the original mesh to the vertices of the cut mesh.
    pub fn lift_vertex_values<T: Copy>(&self, original: &[T]) -> Vec<T> {
        self.vertex_to_original
            .iter()
            .map(|&v| original[v as usize])
            .collect()
    }
}

/// The Euler characteristic of a mesh after cutting, from the Euler characteristic before:
/// `χ' = χ + a − c`, where `a` is the number of vertex copies the cut adds (`Σ_v (s_v − 1)` with
/// `s_v` the number of sectors at `v`) and `c` the number of interior edges that are cut. Faces
/// are unchanged, every cut edge becomes two edges.
///
/// [formula: cut-euler-characteristic]
#[inline]
pub fn cut_euler_characteristic(
    euler_characteristic: i64,
    added_vertices: usize,
    cut_edges: usize,
) -> i64 {
    euler_characteristic + added_vertices as i64 - cut_edges as i64
}

/// Cuts `mesh` along the given edges (edge ids as in [`TriMesh::edges`]).
///
/// Repeated ids are allowed and count once. Boundary edges in the list are ignored and counted
/// in [`CutMesh::ignored_boundary_edges`]. The result is always a valid manifold mesh: cutting
/// cannot create non-manifold vertices, because every sector becomes a vertex of its own.
///
/// One kind of cut has no representation as an indexed mesh and is rejected: a single edge
/// whose two end points are interior vertices without any other cut edge. Neither end point is
/// split by such a slit, so its two sides would be the same vertex pair.
///
/// # Errors
/// [`BffError::EdgeOutOfRange`] for an id that is not an edge of the mesh;
/// [`BffError::IsolatedCutEdge`] for an isolated slit.
pub fn cut_mesh(mesh: &TriMesh, edges: &[u32]) -> Result<CutMesh, BffError> {
    let n_edges = mesh.n_edges();
    let mut is_cut = vec![false; n_edges];
    let mut on_boundary_requested = vec![false; 0];
    let mut ignored_boundary_edges = 0usize;
    for (position, &e) in edges.iter().enumerate() {
        if e as usize >= n_edges {
            return Err(BffError::EdgeOutOfRange {
                what: "cut edges",
                position,
                edge: e,
                n_edges,
            });
        }
        if mesh.is_boundary_edge(e as usize) {
            if on_boundary_requested.is_empty() {
                on_boundary_requested = vec![false; n_edges];
            }
            if !on_boundary_requested[e as usize] {
                on_boundary_requested[e as usize] = true;
                ignored_boundary_edges += 1;
            }
        } else {
            is_cut[e as usize] = true;
        }
    }
    cut_mesh_flags(mesh, &is_cut, ignored_boundary_edges)
}

/// [`cut_mesh`] for a per-edge flag array in which only interior edges are set.
pub(crate) fn cut_mesh_flags(
    mesh: &TriMesh,
    is_cut: &[bool],
    ignored_boundary_edges: usize,
) -> Result<CutMesh, BffError> {
    cut_mesh_impl(mesh, is_cut, ignored_boundary_edges, None)
}

/// [`cut_mesh_flags`] that tolerates isolated slits: they are left uncut and listed in
/// `slits`. The topology of the result then differs from the true one by one edge and one
/// boundary loop per slit, in the component that contains it.
pub(crate) fn cut_mesh_lenient(
    mesh: &TriMesh,
    is_cut: &[bool],
    slits: &mut Vec<u32>,
) -> Result<CutMesh, BffError> {
    cut_mesh_impl(mesh, is_cut, 0, Some(slits))
}

fn cut_mesh_impl(
    mesh: &TriMesh,
    is_cut: &[bool],
    ignored_boundary_edges: usize,
    mut slits: Option<&mut Vec<u32>>,
) -> Result<CutMesh, BffError> {
    let n = mesh.n_vertices();
    let mut cut_edges: Vec<u32> = (0..mesh.n_edges() as u32)
        .filter(|&e| is_cut[e as usize])
        .collect();

    // The vertex of the cut mesh at every corner.
    let mut corner_vertex = vec![NONE; mesh.n_halfedges()];
    let mut vertex_to_original: Vec<u32> = (0..n as u32).collect();
    let mut fan: Vec<usize> = Vec::new();
    for v in 0..n {
        fan.clear();
        fan.extend(mesh.outgoing_halfedges(v));
        // Walking counter-clockwise from fan[i - 1] to fan[i] crosses the edge of fan[i]; a new
        // sector starts there if that edge is cut. At an interior vertex the walk is cyclic:
        // from the last halfedge back to fan[0] it crosses the edge of fan[0].
        let wraps = !mesh.is_boundary_vertex(v) && !is_cut[mesh.edge(fan[0])];
        let mut sector = 0usize;
        let first_new = vertex_to_original.len();
        for (i, &h) in fan.iter().enumerate() {
            if i > 0 && is_cut[mesh.edge(h)] {
                sector += 1;
                vertex_to_original.push(v as u32);
            }
            corner_vertex[h] = if sector == 0 {
                v as u32
            } else {
                (first_new + sector - 1) as u32
            };
        }
        if wraps && sector > 0 {
            // The last sector continues across fan[0] into the first one: merge them.
            let last = (first_new + sector - 1) as u32;
            for &h in fan.iter().rev() {
                if corner_vertex[h] != last {
                    break;
                }
                corner_vertex[h] = v as u32;
            }
            vertex_to_original.pop();
        }
    }

    // Both sides of a cut edge must differ in at least one end point.
    let mut found_slit = false;
    for &e in &cut_edges {
        let h = mesh.edge_halfedge(e as usize);
        let t = mesh.twin(h).expect("only interior edges are flagged");
        if corner_vertex[h] == corner_vertex[mesh.next(t)]
            && corner_vertex[mesh.next(h)] == corner_vertex[t]
        {
            match slits.as_deref_mut() {
                Some(list) => {
                    list.push(e);
                    found_slit = true;
                }
                None => {
                    return Err(BffError::IsolatedCutEdge {
                        edge: e,
                        vertices: mesh.edges()[e as usize],
                    });
                }
            }
        }
    }
    if found_slit && let Some(list) = slits.as_deref() {
        cut_edges.retain(|e| !list.contains(e));
    }

    let faces: Vec<[u32; 3]> = (0..mesh.n_faces())
        .map(|f| {
            [
                corner_vertex[3 * f],
                corner_vertex[3 * f + 1],
                corner_vertex[3 * f + 2],
            ]
        })
        .collect();
    let cut = TriMesh::new(vertex_to_original.len(), faces)?;
    let edge_to_original: Vec<u32> = (0..cut.n_edges())
        .map(|e| mesh.edge(cut.edge_halfedge(e)) as u32)
        .collect();
    Ok(CutMesh {
        mesh: cut,
        vertex_to_original,
        edge_to_original,
        cut_edges,
        ignored_boundary_edges,
        n_original_vertices: n,
    })
}

/// The interior edges across which a per-corner attribute is discontinuous: edge `ab` with
/// faces `f` and `g` is listed if the value at the corner of `f` at `a` differs from the value
/// at the corner of `g` at `a`, or likewise at `b`. `corner_values[3 f + k]` belongs to the
/// corner of face `f` at `faces[f][k]`.
///
/// Applied to per-corner UV coordinates this recovers the seams of an existing UV map. A cut
/// edge both of whose end points carry a single value on both sides (an isolated slit) cannot
/// be detected this way.
///
/// # Errors
/// [`BffError::LengthMismatch`] if there is not one value per corner.
pub fn discontinuous_edges<T: PartialEq>(
    mesh: &TriMesh,
    corner_values: &[T],
) -> Result<Vec<u32>, BffError> {
    if corner_values.len() != mesh.n_halfedges() {
        return Err(BffError::LengthMismatch {
            what: "corner values",
            expected: mesh.n_halfedges(),
            got: corner_values.len(),
        });
    }
    let mut edges = Vec::new();
    for e in 0..mesh.n_edges() {
        let h = mesh.edge_halfedge(e);
        let Some(t) = mesh.twin(h) else {
            continue;
        };
        // h runs a -> b in one face, t runs b -> a in the other.
        if corner_values[h] != corner_values[mesh.next(t)]
            || corner_values[mesh.next(h)] != corner_values[t]
        {
            edges.push(e as u32);
        }
    }
    Ok(edges)
}

/// Topology of one connected component of a mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComponentTopology {
    /// Lowest-numbered face of the component.
    pub first_face: u32,
    /// Number of vertices.
    pub n_vertices: usize,
    /// Number of edges.
    pub n_edges: usize,
    /// Number of faces.
    pub n_faces: usize,
    /// Number of boundary loops.
    pub boundary_loops: usize,
    /// Euler characteristic `V − E + F`.
    pub euler_characteristic: i64,
    /// Genus `(2 − b − χ) / 2`.
    pub genus: usize,
}

impl ComponentTopology {
    /// Whether the component is a topological disk: genus 0 with exactly one boundary loop,
    /// equivalently (for a connected orientable surface) Euler characteristic 1.
    ///
    /// [formula: genus-from-euler]
    #[inline]
    pub fn is_disk(&self) -> bool {
        self.boundary_loops == 1 && self.euler_characteristic == 1
    }

    /// Whether the component has no boundary.
    #[inline]
    pub fn is_closed(&self) -> bool {
        self.boundary_loops == 0
    }
}

/// The topology of every connected component of `mesh`, in component order (see
/// [`TriMesh::face_components`]).
///
/// [formula: euler-characteristic]
/// [formula: genus-from-euler]
pub fn component_topology(mesh: &TriMesh) -> Vec<ComponentTopology> {
    let face_component = mesh.face_components();
    let mut out = vec![
        ComponentTopology {
            first_face: NONE,
            n_vertices: 0,
            n_edges: 0,
            n_faces: 0,
            boundary_loops: 0,
            euler_characteristic: 0,
            genus: 0,
        };
        mesh.n_components()
    ];
    for (f, &c) in face_component.iter().enumerate() {
        let t = &mut out[c as usize];
        if t.first_face == NONE {
            t.first_face = f as u32;
        }
        t.n_faces += 1;
    }
    for v in 0..mesh.n_vertices() {
        out[face_component[mesh.outgoing_halfedge(v) / 3] as usize].n_vertices += 1;
    }
    for e in 0..mesh.n_edges() {
        out[face_component[mesh.edge_halfedge(e) / 3] as usize].n_edges += 1;
    }
    for boundary in mesh.boundary_loops() {
        let h = mesh.outgoing_halfedge(boundary[0] as usize);
        out[face_component[h / 3] as usize].boundary_loops += 1;
    }
    for t in &mut out {
        t.euler_characteristic = t.n_vertices as i64 - t.n_edges as i64 + t.n_faces as i64;
        let twice = 2 - t.boundary_loops as i64 - t.euler_characteristic;
        t.genus = (twice.max(0) / 2) as usize;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators;

    #[test]
    fn cutting_nothing_is_the_identity() {
        let g = generators::icosphere(1).unwrap();
        let mesh = g.topology().unwrap();
        let cut = cut_mesh(&mesh, &[]).unwrap();
        assert_eq!(cut.mesh, mesh);
        assert_eq!(cut.n_added_vertices(), 0);
        assert!(cut.cut_edges.is_empty());
        assert_eq!(
            cut.edge_to_original,
            (0..mesh.n_edges() as u32).collect::<Vec<_>>()
        );
    }

    #[test]
    fn a_two_edge_slit_in_a_sphere_makes_a_disk() {
        let g = generators::icosphere(0).unwrap();
        let mesh = g.topology().unwrap();
        // Two edges sharing vertex 0.
        let mut at_zero = mesh.outgoing_halfedges(0).map(|h| mesh.edge(h) as u32);
        let (e0, e1) = (at_zero.next().unwrap(), at_zero.next().unwrap());
        let cut = cut_mesh(&mesh, &[e0, e1]).unwrap();
        // Only the middle vertex is duplicated; two edges become four.
        assert_eq!(cut.n_added_vertices(), 1);
        assert_eq!(cut.vertex_to_original[mesh.n_vertices()], 0);
        assert_eq!(cut.mesh.n_edges(), mesh.n_edges() + 2);
        assert_eq!(cut.mesh.boundary_loops().len(), 1);
        assert_eq!(cut.mesh.boundary_loops()[0].len(), 4);
        let t = component_topology(&cut.mesh);
        assert_eq!(t.len(), 1);
        assert!(t[0].is_disk());
        assert_eq!(
            cut.mesh.euler_characteristic(),
            cut_euler_characteristic(mesh.euler_characteristic(), 1, 2)
        );
        // A single edge is a slit whose sides cannot be told apart.
        let err = cut_mesh(&mesh, &[e0]).unwrap_err();
        assert!(matches!(err, BffError::IsolatedCutEdge { edge, .. } if edge == e0), "{err}");
    }

    #[test]
    fn cutting_every_edge_isolates_every_face() {
        let g = generators::torus(5, 4, 2.0, 0.5).unwrap();
        let mesh = g.topology().unwrap();
        let all: Vec<u32> = (0..mesh.n_edges() as u32).collect();
        let cut = cut_mesh(&mesh, &all).unwrap();
        assert_eq!(cut.mesh.n_components(), mesh.n_faces());
        assert_eq!(cut.mesh.n_vertices(), 3 * mesh.n_faces());
        assert!(component_topology(&cut.mesh).iter().all(|t| t.is_disk()));
        assert_eq!(
            cut.mesh.euler_characteristic(),
            cut_euler_characteristic(
                mesh.euler_characteristic(),
                cut.n_added_vertices(),
                cut.cut_edges.len()
            )
        );
        for (v, &o) in cut.vertex_to_original.iter().enumerate() {
            for h in cut.mesh.outgoing_halfedges(v) {
                assert_eq!(mesh.from_vertex(h), o as usize);
            }
        }
    }

    #[test]
    fn boundary_edges_are_ignored_and_bad_ids_rejected() {
        let g = generators::grid_disk(3, 3).unwrap();
        let mesh = g.topology().unwrap();
        let boundary_edge = (0..mesh.n_edges())
            .find(|&e| mesh.is_boundary_edge(e))
            .unwrap() as u32;
        let cut = cut_mesh(&mesh, &[boundary_edge, boundary_edge]).unwrap();
        assert_eq!(cut.ignored_boundary_edges, 1);
        assert_eq!(cut.mesh, mesh);
        let err = cut_mesh(&mesh, &[mesh.n_edges() as u32]).unwrap_err();
        assert!(matches!(err, BffError::EdgeOutOfRange { .. }), "{err}");
    }

    #[test]
    fn cuts_are_recovered_from_per_corner_values() {
        let g = generators::icosphere(2).unwrap();
        let mesh = g.topology().unwrap();
        // A path of three edges starting at vertex 0.
        let mut path = Vec::new();
        let mut v = 0usize;
        for _ in 0..3 {
            let h = mesh
                .outgoing_halfedges(v)
                .find(|&h| !path.contains(&(mesh.edge(h) as u32)))
                .unwrap();
            path.push(mesh.edge(h) as u32);
            v = mesh.to_vertex(h);
        }
        path.sort_unstable();
        let cut = cut_mesh(&mesh, &path).unwrap();
        // The vertex of the cut mesh at every corner is a per-corner value that is
        // discontinuous exactly across the cut.
        let corner_vertex: Vec<u32> = cut.mesh.faces().iter().flatten().copied().collect();
        assert_eq!(discontinuous_edges(&mesh, &corner_vertex).unwrap(), path);
        assert!(discontinuous_edges(&mesh, &corner_vertex[1..]).is_err());
    }

    #[test]
    fn a_chord_splits_a_disk_in_two() {
        // 2 x 1 grid: cutting the middle vertical edge gives two squares.
        let g = generators::grid_disk(2, 1).unwrap();
        let mesh = g.topology().unwrap();
        // Vertex (i, j) has index i * 2 + j; the middle edge joins (1, 0) and (1, 1).
        let e = mesh.find_edge(2, 3).unwrap() as u32;
        let cut = cut_mesh(&mesh, &[e]).unwrap();
        assert_eq!(cut.mesh.n_components(), 2);
        assert_eq!(cut.n_added_vertices(), 2);
        assert_eq!(cut.copies_per_vertex()[2], 2);
        assert_eq!(cut.copies_per_vertex()[3], 2);
        assert!(component_topology(&cut.mesh).iter().all(|t| t.is_disk()));
        // Both sides of the cut map to the same original edge.
        let sides = cut
            .edge_to_original
            .iter()
            .filter(|&&o| o == e)
            .count();
        assert_eq!(sides, 2);
    }
}
