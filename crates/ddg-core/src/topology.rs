//! Validated, oriented, manifold triangle mesh connectivity.
//!
//! # Indexing conventions (the contract for every other crate)
//!
//! * Vertices are `0..n_vertices`, faces are `0..n_faces`.
//! * **Halfedge** `h = 3 f + k` (`k ∈ {0,1,2}`) lives in face `f` and runs from
//!   `faces[f][k]` to `faces[f][(k + 1) % 3]`. There are exactly `3 n_faces` halfedges;
//!   boundary edges have one halfedge only (there are no "ghost" halfedges outside the mesh).
//! * **Corner** `c = 3 f + k` is the corner of face `f` at vertex `faces[f][k]`, so a corner
//!   and the halfedge leaving it share one id. The halfedge *opposite* corner `c` is
//!   [`TriMesh::next`]`(c)`; the corner opposite halfedge `h` is [`TriMesh::prev`]`(h)`.
//! * **Edges** are the unique unordered vertex pairs, stored as `[min, max]` and numbered in
//!   lexicographic order of that pair. The numbering depends only on the edge set.
//! * Faces are counter-clockwise when seen from the side the surface normal points to. Walking
//!   around a vertex with [`TriMesh::outgoing_halfedges`] is counter-clockwise in that sense, and
//!   a boundary loop is traversed with the surface on its left.
//!
//! Ids are passed to and returned from methods as `usize`; stored arrays use `u32`.
//!
//! # What a `TriMesh` guarantees
//!
//! [`TriMesh::new`] accepts a face list only if
//!
//! 1. every index is in range and no face repeats a vertex,
//! 2. no two faces have the same three vertices,
//! 3. every edge has one or two incident faces,
//! 4. the two faces of an interior edge traverse it in opposite directions (consistent orientation),
//! 5. the faces around every vertex form exactly one fan (an open fan at a boundary vertex, a
//!    closed one at an interior vertex), and
//! 6. every vertex is used by at least one face.
//!
//! Point 6 is a deliberate policy: an unreferenced vertex would give the cotan Laplacian an empty
//! row. Use [`compact_vertices`] to drop such vertices and keep an index map. Inconsistently
//! oriented input can be repaired with [`orient_faces`]. Neither repair happens implicitly.

use std::fmt;

/// Sentinel used in stored `u32` index arrays for "no element".
pub const NONE: u32 = u32::MAX;

/// Why a face list is not a valid oriented manifold triangle mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshError {
    /// The face list is empty.
    NoFaces,
    /// More vertices or faces than `u32` ids can address.
    TooLarge {
        /// Number of vertices requested.
        n_vertices: usize,
        /// Number of faces requested.
        n_faces: usize,
    },
    /// A face references a vertex index `>= n_vertices`.
    VertexIndexOutOfRange {
        /// Face containing the bad index.
        face: u32,
        /// Corner slot (0, 1 or 2) of the bad index.
        corner: u8,
        /// The offending index.
        index: u32,
        /// Number of vertices in the mesh.
        n_vertices: usize,
    },
    /// A face uses the same vertex more than once.
    RepeatedVertexInFace {
        /// The offending face.
        face: u32,
        /// Its three vertex indices.
        vertices: [u32; 3],
    },
    /// Two faces have the same three vertices (in any order or orientation).
    DuplicateFace {
        /// The lower-numbered face.
        first: u32,
        /// The higher-numbered face.
        second: u32,
        /// The shared vertices, sorted.
        vertices: [u32; 3],
    },
    /// An edge has more than two incident faces.
    NonManifoldEdge {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// All faces containing the edge, ascending.
        faces: Vec<u32>,
    },
    /// The two faces of an edge traverse it in the same direction.
    InconsistentOrientation {
        /// The edge as `[min, max]`.
        edge: [u32; 2],
        /// The two faces, ascending.
        faces: [u32; 2],
    },
    /// The faces around a vertex form more than one fan.
    NonManifoldVertex {
        /// The offending vertex.
        vertex: u32,
        /// Number of separate fans of faces around it.
        fans: u32,
    },
    /// A vertex is used by no face.
    UnreferencedVertex {
        /// The lowest-numbered unreferenced vertex.
        vertex: u32,
        /// Total number of unreferenced vertices.
        count: usize,
    },
    /// A connected component admits no consistent orientation (reported by [`orient_faces`]).
    NonOrientable {
        /// Index of the component (components are numbered by their lowest face).
        component: u32,
        /// The edge across which the conflict was detected, as `[min, max]`.
        edge: [u32; 2],
        /// The two faces that cannot both be oriented consistently.
        faces: [u32; 2],
    },
    /// A face index passed to [`TriMesh::submesh`] is out of range or listed twice.
    BadFaceSelection {
        /// Position in the selection list.
        position: usize,
        /// The offending face index.
        face: u32,
        /// Number of faces in the mesh.
        n_faces: usize,
    },
}

impl fmt::Display for MeshError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MeshError::NoFaces => write!(f, "the mesh has no faces"),
            MeshError::TooLarge {
                n_vertices,
                n_faces,
            } => write!(
                f,
                "mesh too large for 32-bit ids: {n_vertices} vertices, {n_faces} faces \
                 (limits: {} vertices, {} faces)",
                NONE - 1,
                (NONE - 1) / 3
            ),
            MeshError::VertexIndexOutOfRange {
                face,
                corner,
                index,
                n_vertices,
            } => {
                if *n_vertices == 0 {
                    write!(
                        f,
                        "face {face} corner {corner} references vertex {index}, but the mesh has no vertices"
                    )
                } else {
                    write!(
                        f,
                        "face {face} corner {corner} references vertex {index}, but the mesh has \
                         {n_vertices} vertices (valid indices 0..={})",
                        n_vertices - 1
                    )
                }
            }
            MeshError::RepeatedVertexInFace {
                face,
                vertices: [a, b, c],
            } => {
                let repeated = if a == b || a == c { a } else { b };
                write!(f, "face {face} = ({a}, {b}, {c}) repeats vertex {repeated}")
            }
            MeshError::DuplicateFace {
                first,
                second,
                vertices: [a, b, c],
            } => write!(
                f,
                "faces {first} and {second} are duplicates: both have the vertices ({a}, {b}, {c})"
            ),
            MeshError::NonManifoldEdge {
                edge: [a, b],
                faces,
            } => {
                write!(f, "edge ({a}, {b}) is shared by {} faces: ", faces.len())?;
                write_list(f, faces)
            }
            MeshError::InconsistentOrientation {
                edge: [a, b],
                faces: [f0, f1],
            } => write!(
                f,
                "faces {f0} and {f1} traverse their shared edge ({a}, {b}) in the same direction: \
                 the mesh is not consistently oriented (orient_faces can repair this)"
            ),
            MeshError::NonManifoldVertex { vertex, fans } => write!(
                f,
                "vertex {vertex} is non-manifold: its incident faces form {fans} separate fans"
            ),
            MeshError::UnreferencedVertex { vertex, count } => write!(
                f,
                "vertex {vertex} is not referenced by any face ({count} unreferenced \
                 {} in total; compact_vertices removes them)",
                if *count == 1 { "vertex" } else { "vertices" }
            ),
            MeshError::NonOrientable {
                component,
                edge: [a, b],
                faces: [f0, f1],
            } => write!(
                f,
                "connected component {component} is not orientable: faces {f0} and {f1} cannot be \
                 oriented consistently across edge ({a}, {b})"
            ),
            MeshError::BadFaceSelection {
                position,
                face,
                n_faces,
            } => {
                if (*face as usize) < *n_faces {
                    write!(
                        f,
                        "face selection entry {position} lists face {face} a second time"
                    )
                } else {
                    write!(
                        f,
                        "face selection entry {position} is face {face}, but the mesh has {n_faces} faces"
                    )
                }
            }
        }
    }
}

impl std::error::Error for MeshError {}

/// Writes at most 16 items of a list, then an ellipsis with the remaining count.
pub(crate) fn write_list<T: fmt::Display>(f: &mut fmt::Formatter<'_>, items: &[T]) -> fmt::Result {
    const MAX: usize = 16;
    for (i, item) in items.iter().take(MAX).enumerate() {
        if i > 0 {
            write!(f, ", ")?;
        }
        write!(f, "{item}")?;
    }
    if items.len() > MAX {
        write!(f, ", ... ({} more)", items.len() - MAX)?;
    }
    Ok(())
}

/// A validated, oriented, manifold triangle mesh (connectivity only, no geometry).
///
/// See the [module documentation](self) for the indexing conventions and guarantees.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TriMesh {
    n_vertices: usize,
    faces: Vec<[u32; 3]>,
    /// Per halfedge: the oppositely directed halfedge in the adjacent face, or `NONE`.
    twin: Vec<u32>,
    /// Per vertex: an outgoing halfedge; the one without a twin if the vertex is on the boundary.
    vertex_halfedge: Vec<u32>,
    /// Per halfedge: its edge.
    halfedge_edge: Vec<u32>,
    /// Per edge: `[min, max]` vertices; edges are sorted lexicographically.
    edges: Vec<[u32; 2]>,
    /// Per edge: the lower-numbered of its halfedges.
    edge_halfedge: Vec<u32>,
    boundary_loops: Vec<Vec<u32>>,
    n_boundary_edges: usize,
    face_component: Vec<u32>,
    n_components: usize,
}

/// A part of a [`TriMesh`] as a mesh of its own, with maps back to the parent.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubMesh {
    /// The extracted mesh, with vertices and faces renumbered from 0.
    pub mesh: TriMesh,
    /// For each vertex of `mesh`, its index in the parent mesh (strictly increasing).
    pub vertex_to_parent: Vec<u32>,
    /// For each face of `mesh`, its index in the parent mesh. Corner slots are unchanged:
    /// corner `k` of face `f` is corner `k` of parent face `face_to_parent[f]`.
    pub face_to_parent: Vec<u32>,
}

/// Result of [`compact_vertices`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompactedVertices {
    /// The faces with renumbered vertex indices (same order, same corner slots).
    pub faces: Vec<[u32; 3]>,
    /// For each kept vertex, its original index (strictly increasing). Its length is the new
    /// vertex count.
    pub new_to_old: Vec<u32>,
    /// For each original vertex, its new index, or [`NONE`] if it was unreferenced and removed.
    pub old_to_new: Vec<u32>,
}

/// Result of [`orient_faces`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrientedFaces {
    /// The consistently oriented faces. A flipped face `[a, b, c]` becomes `[a, c, b]`:
    /// corner slot 0 is unchanged and slots 1 and 2 are swapped.
    pub faces: Vec<[u32; 3]>,
    /// For each face, whether it was flipped.
    pub flipped: Vec<bool>,
    /// Number of flipped faces.
    pub n_flipped: usize,
    /// Number of connected components (faces connected across shared edges).
    pub n_components: usize,
}

#[inline]
fn edge_key(a: u32, b: u32) -> u64 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    (u64::from(lo) << 32) | u64::from(hi)
}

#[inline]
fn key_edge(key: u64) -> [u32; 2] {
    [(key >> 32) as u32, key as u32]
}

/// Checks size limits, index ranges and repeated vertices.
fn check_faces_basic(n_vertices: usize, faces: &[[u32; 3]]) -> Result<(), MeshError> {
    if n_vertices >= NONE as usize || faces.len() > ((NONE - 1) / 3) as usize {
        return Err(MeshError::TooLarge {
            n_vertices,
            n_faces: faces.len(),
        });
    }
    for (f, face) in faces.iter().enumerate() {
        for (k, &v) in face.iter().enumerate() {
            if v as usize >= n_vertices {
                return Err(MeshError::VertexIndexOutOfRange {
                    face: f as u32,
                    corner: k as u8,
                    index: v,
                    n_vertices,
                });
            }
        }
        if face[0] == face[1] || face[1] == face[2] || face[0] == face[2] {
            return Err(MeshError::RepeatedVertexInFace {
                face: f as u32,
                vertices: *face,
            });
        }
    }
    Ok(())
}

/// All halfedges as `(undirected edge key, halfedge id)`, sorted. Halfedges of one edge are
/// adjacent and in ascending halfedge order.
fn sorted_halfedge_keys(faces: &[[u32; 3]]) -> Vec<(u64, u32)> {
    let mut keys = Vec::with_capacity(faces.len() * 3);
    for (f, face) in faces.iter().enumerate() {
        for k in 0..3 {
            keys.push((edge_key(face[k], face[(k + 1) % 3]), (3 * f + k) as u32));
        }
    }
    keys.sort_unstable();
    keys
}

/// Returns the first edge with more than two incident faces as an error, if any.
fn check_edge_manifold(keys: &[(u64, u32)]) -> Result<(), MeshError> {
    let mut i = 0;
    while i < keys.len() {
        let mut j = i + 1;
        while j < keys.len() && keys[j].0 == keys[i].0 {
            j += 1;
        }
        if j - i > 2 {
            return Err(MeshError::NonManifoldEdge {
                edge: key_edge(keys[i].0),
                faces: keys[i..j].iter().map(|&(_, h)| h / 3).collect(),
            });
        }
        i = j;
    }
    Ok(())
}

/// Removes vertices that no face references and renumbers the rest, keeping their order.
///
/// This is the explicit companion to the "no unreferenced vertices" rule of [`TriMesh::new`].
/// Per-vertex data of the original mesh maps to the compacted mesh through
/// [`CompactedVertices::new_to_old`].
///
/// # Errors
/// [`MeshError::VertexIndexOutOfRange`] if a face references an index `>= n_vertices`, and
/// [`MeshError::TooLarge`] if the counts exceed the `u32` id range.
pub fn compact_vertices(
    n_vertices: usize,
    faces: &[[u32; 3]],
) -> Result<CompactedVertices, MeshError> {
    if n_vertices >= NONE as usize {
        return Err(MeshError::TooLarge {
            n_vertices,
            n_faces: faces.len(),
        });
    }
    let mut old_to_new = vec![NONE; n_vertices];
    for (f, face) in faces.iter().enumerate() {
        for (k, &v) in face.iter().enumerate() {
            if v as usize >= n_vertices {
                return Err(MeshError::VertexIndexOutOfRange {
                    face: f as u32,
                    corner: k as u8,
                    index: v,
                    n_vertices,
                });
            }
            old_to_new[v as usize] = 0;
        }
    }
    let mut new_to_old = Vec::new();
    for (old, slot) in old_to_new.iter_mut().enumerate() {
        if *slot != NONE {
            *slot = new_to_old.len() as u32;
            new_to_old.push(old as u32);
        }
    }
    let faces = faces
        .iter()
        .map(|face| face.map(|v| old_to_new[v as usize]))
        .collect();
    Ok(CompactedVertices {
        faces,
        new_to_old,
        old_to_new,
    })
}

/// Re-orients faces so that every connected component is consistently oriented.
///
/// Faces are connected when they share an edge. Within each component the orientation that
/// flips the fewest faces is chosen; on a tie the component's lowest-numbered face keeps its
/// orientation. Which side is "outside" is not decided here: that needs geometry.
///
/// A flipped face `[a, b, c]` becomes `[a, c, b]`, so per-corner data of a flipped face must
/// have its slots 1 and 2 swapped by the caller (see [`OrientedFaces::flipped`]).
///
/// The result is not otherwise validated; pass it to [`TriMesh::new`].
///
/// # Errors
/// * [`MeshError::VertexIndexOutOfRange`], [`MeshError::RepeatedVertexInFace`],
///   [`MeshError::TooLarge`] for malformed faces;
/// * [`MeshError::NonManifoldEdge`] if an edge has more than two faces (orientation across such
///   an edge is ambiguous);
/// * [`MeshError::NonOrientable`] if a component (for example a Möbius strip) has no consistent
///   orientation.
pub fn orient_faces(n_vertices: usize, faces: &[[u32; 3]]) -> Result<OrientedFaces, MeshError> {
    check_faces_basic(n_vertices, faces)?;
    let keys = sorted_halfedge_keys(faces);
    check_edge_manifold(&keys)?;

    let n_halfedges = faces.len() * 3;
    // Neighbour across each halfedge, regardless of direction.
    let mut mate = vec![NONE; n_halfedges];
    let mut i = 0;
    while i + 1 < keys.len() {
        if keys[i].0 == keys[i + 1].0 {
            mate[keys[i].1 as usize] = keys[i + 1].1;
            mate[keys[i + 1].1 as usize] = keys[i].1;
            i += 2;
        } else {
            i += 1;
        }
    }
    let origin = |h: u32| faces[(h / 3) as usize][(h % 3) as usize];

    const UNSEEN: u8 = 2;
    let mut state = vec![UNSEEN; faces.len()]; // 0 = keep, 1 = flip
    let mut n_components = 0usize;
    let mut stack: Vec<u32> = Vec::new();
    let mut members: Vec<u32> = Vec::new();
    for seed in 0..faces.len() {
        if state[seed] != UNSEEN {
            continue;
        }
        state[seed] = 0;
        stack.push(seed as u32);
        members.clear();
        while let Some(f) = stack.pop() {
            members.push(f);
            for k in 0..3 {
                let h = 3 * f + k;
                let m = mate[h as usize];
                if m == NONE {
                    continue;
                }
                let g = m / 3;
                // As stored, the two faces agree across the edge iff they traverse it in
                // opposite directions.
                let same_direction = origin(h) == origin(m);
                let wanted = state[f as usize] ^ u8::from(same_direction);
                if state[g as usize] == UNSEEN {
                    state[g as usize] = wanted;
                    stack.push(g);
                } else if state[g as usize] != wanted {
                    let a = origin(h);
                    let b = faces[f as usize][((k + 1) % 3) as usize];
                    return Err(MeshError::NonOrientable {
                        component: n_components as u32,
                        edge: key_edge(edge_key(a, b)),
                        faces: [f.min(g), f.max(g)],
                    });
                }
            }
        }
        let flips = members.iter().filter(|&&f| state[f as usize] == 1).count();
        if 2 * flips > members.len() {
            for &f in &members {
                state[f as usize] ^= 1;
            }
        }
        n_components += 1;
    }

    let flipped: Vec<bool> = state.iter().map(|&s| s == 1).collect();
    let n_flipped = flipped.iter().filter(|&&b| b).count();
    let faces = faces
        .iter()
        .zip(&flipped)
        .map(|(&[a, b, c], &flip)| if flip { [a, c, b] } else { [a, b, c] })
        .collect();
    Ok(OrientedFaces {
        faces,
        flipped,
        n_flipped,
        n_components,
    })
}

impl TriMesh {
    /// Builds and validates a mesh from `n_vertices` and a triangle list.
    ///
    /// The faces are stored exactly as given (no reordering, no re-orientation).
    ///
    /// # Errors
    /// Checks run in this order and the first failure is returned:
    /// [`MeshError::NoFaces`], [`MeshError::TooLarge`], [`MeshError::VertexIndexOutOfRange`],
    /// [`MeshError::RepeatedVertexInFace`], [`MeshError::DuplicateFace`],
    /// [`MeshError::NonManifoldEdge`], [`MeshError::InconsistentOrientation`], then per vertex in
    /// ascending order [`MeshError::NonManifoldVertex`] or [`MeshError::UnreferencedVertex`].
    pub fn new(n_vertices: usize, faces: Vec<[u32; 3]>) -> Result<TriMesh, MeshError> {
        if faces.is_empty() {
            return Err(MeshError::NoFaces);
        }
        check_faces_basic(n_vertices, &faces)?;

        // Duplicate faces (same vertex set, any orientation).
        {
            let mut sorted: Vec<([u32; 3], u32)> = faces
                .iter()
                .enumerate()
                .map(|(f, face)| {
                    let mut s = *face;
                    s.sort_unstable();
                    (s, f as u32)
                })
                .collect();
            sorted.sort_unstable();
            if let Some(w) = sorted.windows(2).find(|w| w[0].0 == w[1].0) {
                return Err(MeshError::DuplicateFace {
                    first: w[0].1,
                    second: w[1].1,
                    vertices: w[0].0,
                });
            }
        }

        let keys = sorted_halfedge_keys(&faces);
        check_edge_manifold(&keys)?;

        let n_halfedges = faces.len() * 3;
        let origin = |h: u32| faces[(h / 3) as usize][(h % 3) as usize];
        let mut twin = vec![NONE; n_halfedges];
        let mut halfedge_edge = vec![NONE; n_halfedges];
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(n_halfedges / 2 + 1);
        let mut edge_halfedge: Vec<u32> = Vec::with_capacity(n_halfedges / 2 + 1);
        let mut n_boundary_edges = 0usize;
        let mut i = 0;
        while i < keys.len() {
            let (key, h0) = keys[i];
            let e = edges.len() as u32;
            edges.push(key_edge(key));
            edge_halfedge.push(h0);
            halfedge_edge[h0 as usize] = e;
            if i + 1 < keys.len() && keys[i + 1].0 == key {
                let h1 = keys[i + 1].1;
                if origin(h0) == origin(h1) {
                    return Err(MeshError::InconsistentOrientation {
                        edge: key_edge(key),
                        faces: [h0 / 3, h1 / 3],
                    });
                }
                twin[h0 as usize] = h1;
                twin[h1 as usize] = h0;
                halfedge_edge[h1 as usize] = e;
                i += 2;
            } else {
                n_boundary_edges += 1;
                i += 1;
            }
        }

        // Fans around vertices. The rotation h -> twin(prev(h)) maps an outgoing halfedge of a
        // vertex to the next outgoing halfedge counter-clockwise. Its orbits are the fans:
        // chains that start at a halfedge without twin (open fans), or cycles (closed fans).
        let prev = |h: u32| if h.is_multiple_of(3) { h + 2 } else { h - 1 };
        let mut fans = vec![0u32; n_vertices];
        let mut vertex_halfedge = vec![NONE; n_vertices];
        let mut visited = vec![false; n_halfedges];
        for h0 in 0..n_halfedges as u32 {
            if twin[h0 as usize] != NONE {
                continue;
            }
            let v = origin(h0) as usize;
            fans[v] += 1;
            if vertex_halfedge[v] == NONE {
                vertex_halfedge[v] = h0;
            }
            let mut h = h0;
            loop {
                visited[h as usize] = true;
                let t = twin[prev(h) as usize];
                if t == NONE {
                    break;
                }
                h = t;
            }
        }
        for h0 in 0..n_halfedges as u32 {
            if visited[h0 as usize] {
                continue;
            }
            let v = origin(h0) as usize;
            fans[v] += 1;
            if vertex_halfedge[v] == NONE {
                vertex_halfedge[v] = h0;
            }
            let mut h = h0;
            loop {
                visited[h as usize] = true;
                let t = twin[prev(h) as usize];
                if t == NONE || t == h0 {
                    break;
                }
                h = t;
            }
        }
        for (v, &count) in fans.iter().enumerate() {
            if count > 1 {
                return Err(MeshError::NonManifoldVertex {
                    vertex: v as u32,
                    fans: count,
                });
            }
            if count == 0 {
                return Err(MeshError::UnreferencedVertex {
                    vertex: v as u32,
                    count: fans.iter().filter(|&&c| c == 0).count(),
                });
            }
        }

        // Boundary loops: from a boundary halfedge a -> b, the next one is the boundary
        // halfedge leaving b, which is vertex_halfedge[b] because b has a single open fan.
        let mut boundary_loops: Vec<Vec<u32>> = Vec::new();
        let mut on_loop = vec![false; n_halfedges];
        for h0 in 0..n_halfedges as u32 {
            if twin[h0 as usize] != NONE || on_loop[h0 as usize] {
                continue;
            }
            let mut vertices = Vec::new();
            let mut h = h0;
            loop {
                on_loop[h as usize] = true;
                vertices.push(origin(h));
                let next = if h % 3 == 2 { h - 2 } else { h + 1 };
                h = vertex_halfedge[origin(next) as usize];
                if h == h0 {
                    break;
                }
            }
            boundary_loops.push(vertices);
        }

        // Connected components of faces, numbered by lowest face.
        let mut face_component = vec![NONE; faces.len()];
        let mut n_components = 0usize;
        let mut stack: Vec<u32> = Vec::new();
        for seed in 0..faces.len() {
            if face_component[seed] != NONE {
                continue;
            }
            face_component[seed] = n_components as u32;
            stack.push(seed as u32);
            while let Some(f) = stack.pop() {
                for k in 0..3 {
                    let t = twin[(3 * f + k) as usize];
                    if t != NONE && face_component[(t / 3) as usize] == NONE {
                        face_component[(t / 3) as usize] = n_components as u32;
                        stack.push(t / 3);
                    }
                }
            }
            n_components += 1;
        }

        Ok(TriMesh {
            n_vertices,
            faces,
            twin,
            vertex_halfedge,
            halfedge_edge,
            edges,
            edge_halfedge,
            boundary_loops,
            n_boundary_edges,
            face_component,
            n_components,
        })
    }

    // ----- counts -----

    /// Number of vertices.
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// Number of faces.
    #[inline]
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of (unique, undirected) edges.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Number of halfedges, which is also the number of corners: `3 * n_faces`.
    #[inline]
    pub fn n_halfedges(&self) -> usize {
        self.faces.len() * 3
    }

    /// Number of boundary edges (edges with a single incident face). Equals the number of
    /// boundary vertices.
    #[inline]
    pub fn n_boundary_edges(&self) -> usize {
        self.n_boundary_edges
    }

    /// Number of connected components.
    #[inline]
    pub fn n_components(&self) -> usize {
        self.n_components
    }

    // ----- raw arrays -----

    /// The faces, exactly as passed to [`TriMesh::new`].
    #[inline]
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// The unique edges as `[min, max]` vertex pairs, sorted lexicographically.
    #[inline]
    pub fn edges(&self) -> &[[u32; 2]] {
        &self.edges
    }

    // ----- halfedges and corners -----

    /// The next halfedge in the same face.
    #[inline]
    pub fn next(&self, h: usize) -> usize {
        if h % 3 == 2 { h - 2 } else { h + 1 }
    }

    /// The previous halfedge in the same face.
    #[inline]
    pub fn prev(&self, h: usize) -> usize {
        if h.is_multiple_of(3) { h + 2 } else { h - 1 }
    }

    /// The oppositely directed halfedge in the neighbouring face, or `None` on the boundary.
    #[inline]
    pub fn twin(&self, h: usize) -> Option<usize> {
        let t = self.twin[h];
        (t != NONE).then_some(t as usize)
    }

    /// The face containing halfedge (or corner) `h`: `h / 3`.
    #[inline]
    pub fn face(&self, h: usize) -> usize {
        h / 3
    }

    /// The vertex halfedge `h` starts at; equivalently the vertex of corner `h`.
    #[inline]
    pub fn from_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][h % 3] as usize
    }

    /// The vertex halfedge `h` points to.
    #[inline]
    pub fn to_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][(h + 1) % 3] as usize
    }

    /// The vertex of the face of `h` that is not on `h`.
    #[inline]
    pub fn opposite_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][(h + 2) % 3] as usize
    }

    /// The edge halfedge `h` belongs to.
    #[inline]
    pub fn edge(&self, h: usize) -> usize {
        self.halfedge_edge[h] as usize
    }

    /// The directed halfedge from `a` to `b`, if it exists. `O(degree(a))`.
    pub fn find_halfedge(&self, a: usize, b: usize) -> Option<usize> {
        if a >= self.n_vertices {
            return None;
        }
        self.outgoing_halfedges(a).find(|&h| self.to_vertex(h) == b)
    }

    // ----- edges -----

    /// The two vertices of edge `e` as `[min, max]`.
    #[inline]
    pub fn edge_vertices(&self, e: usize) -> [usize; 2] {
        let [a, b] = self.edges[e];
        [a as usize, b as usize]
    }

    /// The lower-numbered halfedge of edge `e`. Its twin, if any, is the other halfedge of the
    /// edge. Its direction can be either `min -> max` or `max -> min`.
    #[inline]
    pub fn edge_halfedge(&self, e: usize) -> usize {
        self.edge_halfedge[e] as usize
    }

    /// Whether edge `e` has only one incident face.
    #[inline]
    pub fn is_boundary_edge(&self, e: usize) -> bool {
        self.twin[self.edge_halfedge[e] as usize] == NONE
    }

    /// The edge between vertices `a` and `b`, if they are adjacent. `O(log n_edges)`.
    pub fn find_edge(&self, a: usize, b: usize) -> Option<usize> {
        if a >= self.n_vertices || b >= self.n_vertices {
            return None;
        }
        let key = if a < b {
            [a as u32, b as u32]
        } else {
            [b as u32, a as u32]
        };
        self.edges.binary_search(&key).ok()
    }

    // ----- vertices -----

    /// An outgoing halfedge of vertex `v`: the first one of [`TriMesh::outgoing_halfedges`].
    /// For a boundary vertex it is the boundary halfedge leaving `v`.
    #[inline]
    pub fn outgoing_halfedge(&self, v: usize) -> usize {
        self.vertex_halfedge[v] as usize
    }

    /// Whether vertex `v` lies on the boundary (its fan of faces is open).
    #[inline]
    pub fn is_boundary_vertex(&self, v: usize) -> bool {
        self.twin[self.vertex_halfedge[v] as usize] == NONE
    }

    /// The halfedges leaving `v`, one per incident face, in counter-clockwise order.
    ///
    /// Each item `h` is also the id of the corner of face `h / 3` at `v`. At a boundary vertex
    /// the iteration starts with the boundary halfedge leaving `v`; note that the last boundary
    /// neighbour is then reached only by an incoming halfedge, see [`TriMesh::vertex_neighbors`].
    pub fn outgoing_halfedges(&self, v: usize) -> OutgoingHalfedges<'_> {
        OutgoingHalfedges {
            mesh: self,
            start: self.vertex_halfedge[v],
            current: self.vertex_halfedge[v],
        }
    }

    /// The vertices adjacent to `v` in counter-clockwise order: `k` of them at an interior
    /// vertex with `k` incident faces, `k + 1` at a boundary vertex (from the boundary neighbour
    /// reached by the outgoing boundary halfedge to the one reached by the incoming one).
    pub fn vertex_neighbors(&self, v: usize) -> VertexNeighbors<'_> {
        VertexNeighbors {
            halfedges: self.outgoing_halfedges(v),
            last: None,
        }
    }

    /// Number of faces incident to `v`.
    pub fn vertex_face_count(&self, v: usize) -> usize {
        self.outgoing_halfedges(v).count()
    }

    /// All boundary vertices, ascending.
    pub fn boundary_vertices(&self) -> Vec<u32> {
        (0..self.n_vertices)
            .filter(|&v| self.is_boundary_vertex(v))
            .map(|v| v as u32)
            .collect()
    }

    /// All interior (non-boundary) vertices, ascending.
    pub fn interior_vertices(&self) -> Vec<u32> {
        (0..self.n_vertices)
            .filter(|&v| !self.is_boundary_vertex(v))
            .map(|v| v as u32)
            .collect()
    }

    // ----- boundary -----

    /// The boundary loops. Each loop lists its vertices in order, with the surface on the left;
    /// the boundary halfedge from `loop[i]` to `loop[(i + 1) % len]` is
    /// [`TriMesh::outgoing_halfedge`]`(loop[i])`. Loops are ordered by their lowest halfedge and
    /// each starts at the origin of that halfedge.
    #[inline]
    pub fn boundary_loops(&self) -> &[Vec<u32>] {
        &self.boundary_loops
    }

    // ----- global topology -----

    /// Euler characteristic `χ = V − E + F`.
    ///
    /// [formula: euler-characteristic]
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices as i64 - self.edges.len() as i64 + self.faces.len() as i64
    }

    /// Total genus: the sum of the genera of all connected components,
    /// `g = (2 C − b − χ) / 2` with `C` components and `b` boundary loops.
    ///
    /// For a connected mesh this is the usual `χ = 2 − 2 g − b`.
    ///
    /// [formula: genus-from-euler]
    pub fn genus(&self) -> usize {
        let twice = 2 * self.n_components as i64
            - self.boundary_loops.len() as i64
            - self.euler_characteristic();
        debug_assert!(
            twice >= 0 && twice % 2 == 0,
            "orientable manifold has integer genus"
        );
        (twice / 2) as usize
    }

    /// For each face, the index of its connected component. Components are numbered in order
    /// of their lowest face.
    #[inline]
    pub fn face_components(&self) -> &[u32] {
        &self.face_component
    }

    /// For each vertex, the index of its connected component (numbering as in
    /// [`TriMesh::face_components`]).
    pub fn vertex_components(&self) -> Vec<u32> {
        (0..self.n_vertices)
            .map(|v| self.face_component[self.vertex_halfedge[v] as usize / 3])
            .collect()
    }

    /// Extracts the given faces as a mesh of their own.
    ///
    /// Faces keep the order of `faces` and their corner slots; vertices keep their relative
    /// order. The selection must itself be a valid mesh: removing faces can create a
    /// non-manifold vertex, which is reported as an error (with indices of the *sub*mesh).
    ///
    /// # Errors
    /// [`MeshError::BadFaceSelection`] for an out-of-range or repeated face index, otherwise
    /// whatever [`TriMesh::new`] reports for the selected faces.
    pub fn submesh(&self, faces: &[u32]) -> Result<SubMesh, MeshError> {
        let mut seen = vec![false; self.faces.len()];
        for (position, &f) in faces.iter().enumerate() {
            if f as usize >= self.faces.len() || seen[f as usize] {
                return Err(MeshError::BadFaceSelection {
                    position,
                    face: f,
                    n_faces: self.faces.len(),
                });
            }
            seen[f as usize] = true;
        }
        let selected: Vec<[u32; 3]> = faces.iter().map(|&f| self.faces[f as usize]).collect();
        let compact = compact_vertices(self.n_vertices, &selected)?;
        let mesh = TriMesh::new(compact.new_to_old.len(), compact.faces)?;
        Ok(SubMesh {
            mesh,
            vertex_to_parent: compact.new_to_old,
            face_to_parent: faces.to_vec(),
        })
    }

    /// Splits the mesh into its connected components, in component order. Within each
    /// component faces and vertices keep their relative order.
    pub fn split_components(&self) -> Vec<SubMesh> {
        let mut face_lists: Vec<Vec<u32>> = vec![Vec::new(); self.n_components];
        for (f, &c) in self.face_component.iter().enumerate() {
            face_lists[c as usize].push(f as u32);
        }
        face_lists
            .iter()
            .map(|faces| {
                self.submesh(faces)
                    .expect("a connected component of a valid mesh is a valid mesh")
            })
            .collect()
    }
}

/// Iterator over the halfedges leaving a vertex; see [`TriMesh::outgoing_halfedges`].
#[derive(Debug, Clone)]
pub struct OutgoingHalfedges<'a> {
    mesh: &'a TriMesh,
    start: u32,
    current: u32,
}

impl Iterator for OutgoingHalfedges<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        if self.current == NONE {
            return None;
        }
        let h = self.current as usize;
        let t = self.mesh.twin[self.mesh.prev(h)];
        self.current = if t == self.start { NONE } else { t };
        Some(h)
    }
}

/// Iterator over the neighbours of a vertex; see [`TriMesh::vertex_neighbors`].
#[derive(Debug, Clone)]
pub struct VertexNeighbors<'a> {
    halfedges: OutgoingHalfedges<'a>,
    /// Origin of the incoming boundary halfedge, pending after the last outgoing halfedge.
    last: Option<usize>,
}

impl Iterator for VertexNeighbors<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.halfedges.next() {
            Some(h) => {
                let mesh = self.halfedges.mesh;
                let p = mesh.prev(h);
                self.last = (mesh.twin[p] == NONE).then(|| mesh.from_vertex(p));
                Some(mesh.to_vertex(h))
            }
            None => self.last.take(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tetrahedron() -> Vec<[u32; 3]> {
        vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]]
    }

    #[test]
    fn tetrahedron_is_a_sphere() {
        let m = TriMesh::new(4, tetrahedron()).unwrap();
        assert_eq!((m.n_vertices(), m.n_edges(), m.n_faces()), (4, 6, 4));
        assert_eq!(m.euler_characteristic(), 2);
        assert_eq!(m.genus(), 0);
        assert_eq!(m.n_components(), 1);
        assert!(m.boundary_loops().is_empty());
        assert_eq!(m.n_boundary_edges(), 0);
        for v in 0..4 {
            assert!(!m.is_boundary_vertex(v));
            assert_eq!(m.vertex_face_count(v), 3);
            let mut nb: Vec<usize> = m.vertex_neighbors(v).collect();
            nb.sort_unstable();
            let expect: Vec<usize> = (0..4).filter(|&w| w != v).collect();
            assert_eq!(nb, expect);
        }
        for h in 0..m.n_halfedges() {
            let t = m.twin(h).unwrap();
            assert_eq!(m.twin(t), Some(h));
            assert_eq!(m.from_vertex(h), m.to_vertex(t));
            assert_eq!(m.to_vertex(h), m.from_vertex(t));
            assert_eq!(m.edge(h), m.edge(t));
            assert_eq!(m.next(m.next(m.next(h))), h);
            assert_eq!(m.prev(m.next(h)), h);
            assert_eq!(m.opposite_vertex(h), m.from_vertex(m.prev(h)));
        }
        assert_eq!(m.edges(), &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]);
        assert_eq!(m.find_edge(3, 1), Some(4));
        assert_eq!(m.find_edge(1, 1), None);
        let h = m.find_halfedge(2, 0).unwrap();
        assert_eq!((m.from_vertex(h), m.to_vertex(h)), (2, 0));
    }

    #[test]
    fn single_triangle_boundary() {
        let m = TriMesh::new(3, vec![[0, 1, 2]]).unwrap();
        assert_eq!(m.euler_characteristic(), 1);
        assert_eq!(m.genus(), 0);
        assert_eq!(m.boundary_loops(), &[vec![0, 1, 2]]);
        assert_eq!(m.n_boundary_edges(), 3);
        for v in 0..3 {
            assert!(m.is_boundary_vertex(v));
            assert_eq!(m.outgoing_halfedges(v).collect::<Vec<_>>(), vec![v]);
        }
        assert_eq!(m.vertex_neighbors(0).collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(m.boundary_vertices(), vec![0, 1, 2]);
        assert!(m.interior_vertices().is_empty());
    }

    #[test]
    fn fan_order_is_counter_clockwise() {
        // A hexagon fan around vertex 0, plus a boundary vertex check.
        let faces: Vec<[u32; 3]> = (0..6).map(|i| [0, 1 + i, 1 + (i + 1) % 6]).collect();
        let m = TriMesh::new(7, faces).unwrap();
        assert!(!m.is_boundary_vertex(0));
        assert_eq!(
            m.vertex_neighbors(0).collect::<Vec<_>>(),
            vec![1, 2, 3, 4, 5, 6]
        );
        // Vertex 1 is on the boundary: faces (0,1,2) and (0,6,1). Counter-clockwise around 1
        // starting at the boundary halfedge 1 -> 2.
        assert!(m.is_boundary_vertex(1));
        assert_eq!(m.vertex_neighbors(1).collect::<Vec<_>>(), vec![2, 0, 6]);
        let h = m.outgoing_halfedge(1);
        assert_eq!((m.from_vertex(h), m.to_vertex(h)), (1, 2));
        assert_eq!(m.twin(h), None);
        assert_eq!(m.boundary_loops(), &[vec![1, 2, 3, 4, 5, 6]]);
        assert_eq!(m.euler_characteristic(), 1);
    }

    #[test]
    fn error_no_faces() {
        let e = TriMesh::new(3, vec![]).unwrap_err();
        assert_eq!(e, MeshError::NoFaces);
        assert_eq!(e.to_string(), "the mesh has no faces");
    }

    #[test]
    fn error_index_out_of_range() {
        let e = TriMesh::new(4, vec![[0, 1, 2], [1, 0, 4]]).unwrap_err();
        assert_eq!(
            e,
            MeshError::VertexIndexOutOfRange {
                face: 1,
                corner: 2,
                index: 4,
                n_vertices: 4
            }
        );
        assert_eq!(
            e.to_string(),
            "face 1 corner 2 references vertex 4, but the mesh has 4 vertices (valid indices 0..=3)"
        );
    }

    #[test]
    fn error_repeated_vertex() {
        let e = TriMesh::new(4, vec![[0, 1, 2], [3, 1, 3]]).unwrap_err();
        assert_eq!(
            e,
            MeshError::RepeatedVertexInFace {
                face: 1,
                vertices: [3, 1, 3]
            }
        );
        assert_eq!(e.to_string(), "face 1 = (3, 1, 3) repeats vertex 3");
    }

    #[test]
    fn error_duplicate_face() {
        // Same orientation and opposite orientation are both duplicates.
        for second in [[1, 2, 0], [2, 1, 0]] {
            let e = TriMesh::new(4, vec![[0, 1, 2], [0, 3, 1], second]).unwrap_err();
            assert_eq!(
                e,
                MeshError::DuplicateFace {
                    first: 0,
                    second: 2,
                    vertices: [0, 1, 2]
                }
            );
            assert_eq!(
                e.to_string(),
                "faces 0 and 2 are duplicates: both have the vertices (0, 1, 2)"
            );
        }
    }

    #[test]
    fn error_edge_with_three_faces() {
        let e = TriMesh::new(5, vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]]).unwrap_err();
        assert_eq!(
            e,
            MeshError::NonManifoldEdge {
                edge: [0, 1],
                faces: vec![0, 1, 2]
            }
        );
        assert_eq!(e.to_string(), "edge (0, 1) is shared by 3 faces: 0, 1, 2");
    }

    #[test]
    fn error_inconsistent_orientation() {
        let e = TriMesh::new(4, vec![[0, 1, 2], [0, 1, 3]]).unwrap_err();
        assert_eq!(
            e,
            MeshError::InconsistentOrientation {
                edge: [0, 1],
                faces: [0, 1]
            }
        );
        assert_eq!(
            e.to_string(),
            "faces 0 and 1 traverse their shared edge (0, 1) in the same direction: the mesh is \
             not consistently oriented (orient_faces can repair this)"
        );
    }

    #[test]
    fn error_non_manifold_vertex() {
        // Two triangles touching only at vertex 2 (a bow tie).
        let e = TriMesh::new(5, vec![[0, 1, 2], [2, 3, 4]]).unwrap_err();
        assert_eq!(e, MeshError::NonManifoldVertex { vertex: 2, fans: 2 });
        assert_eq!(
            e.to_string(),
            "vertex 2 is non-manifold: its incident faces form 2 separate fans"
        );

        // Two tetrahedra sharing one vertex: two closed fans.
        let mut faces = tetrahedron();
        faces.extend(
            tetrahedron()
                .iter()
                .map(|f| f.map(|v| if v == 0 { 0 } else { v + 3 })),
        );
        let e = TriMesh::new(7, faces).unwrap_err();
        assert_eq!(e, MeshError::NonManifoldVertex { vertex: 0, fans: 2 });

        // A closed fan and an open fan at the same vertex.
        let mut faces = tetrahedron();
        faces.push([0, 4, 5]);
        let e = TriMesh::new(6, faces).unwrap_err();
        assert_eq!(e, MeshError::NonManifoldVertex { vertex: 0, fans: 2 });
    }

    #[test]
    fn error_unreferenced_vertex_and_compaction() {
        let faces = vec![[0, 2, 5]];
        let e = TriMesh::new(6, faces.clone()).unwrap_err();
        assert_eq!(
            e,
            MeshError::UnreferencedVertex {
                vertex: 1,
                count: 3
            }
        );
        assert_eq!(
            e.to_string(),
            "vertex 1 is not referenced by any face (3 unreferenced vertices in total; \
             compact_vertices removes them)"
        );
        let c = compact_vertices(6, &faces).unwrap();
        assert_eq!(c.faces, vec![[0, 1, 2]]);
        assert_eq!(c.new_to_old, vec![0, 2, 5]);
        assert_eq!(c.old_to_new, vec![0, NONE, 1, NONE, NONE, 2]);
        TriMesh::new(c.new_to_old.len(), c.faces).unwrap();
        let e = compact_vertices(2, &faces).unwrap_err();
        assert_eq!(
            e,
            MeshError::VertexIndexOutOfRange {
                face: 0,
                corner: 1,
                index: 2,
                n_vertices: 2
            }
        );
    }

    #[test]
    fn orient_repairs_flipped_faces() {
        let mut faces = tetrahedron();
        faces[2] = [faces[2][0], faces[2][2], faces[2][1]];
        assert!(matches!(
            TriMesh::new(4, faces.clone()),
            Err(MeshError::InconsistentOrientation { .. })
        ));
        let o = orient_faces(4, &faces).unwrap();
        assert_eq!(o.flipped, vec![false, false, true, false]);
        assert_eq!(o.n_flipped, 1);
        assert_eq!(o.n_components, 1);
        assert_eq!(o.faces, tetrahedron());
        TriMesh::new(4, o.faces).unwrap();
    }

    #[test]
    fn orient_flips_the_minority() {
        // Three of four faces are "wrong" relative to face 0: face 0 is flipped instead.
        let t = tetrahedron();
        let faces = vec![[t[0][0], t[0][2], t[0][1]], t[1], t[2], t[3]];
        let o = orient_faces(4, &faces).unwrap();
        assert_eq!(o.flipped, vec![true, false, false, false]);
        assert_eq!(o.faces, t);
    }

    #[test]
    fn orient_rejects_moebius_strip() {
        // A Möbius strip with 5 triangles on 5 vertices: faces (i, i+1, i+2) mod 5.
        let faces: Vec<[u32; 3]> = (0..5).map(|i| [i, (i + 1) % 5, (i + 2) % 5]).collect();
        let e = orient_faces(5, &faces).unwrap_err();
        match &e {
            MeshError::NonOrientable { component: 0, .. } => {}
            other => panic!("unexpected error {other:?}"),
        }
        let msg = e.to_string();
        assert!(
            msg.starts_with("connected component 0 is not orientable: faces "),
            "{msg}"
        );
        assert!(
            msg.contains("cannot be oriented consistently across edge ("),
            "{msg}"
        );
    }

    #[test]
    fn orient_reports_three_faces_on_an_edge() {
        let e = orient_faces(5, &[[0, 1, 2], [1, 0, 3], [0, 1, 4]]).unwrap_err();
        assert_eq!(
            e,
            MeshError::NonManifoldEdge {
                edge: [0, 1],
                faces: vec![0, 1, 2]
            }
        );
    }

    #[test]
    fn components_and_submeshes() {
        // A tetrahedron (vertices 0..4) interleaved with a single triangle (vertices 4..7).
        let t = tetrahedron();
        let faces = vec![t[0], [4, 5, 6], t[1], t[2], t[3]];
        let m = TriMesh::new(7, faces).unwrap();
        assert_eq!(m.n_components(), 2);
        assert_eq!(m.face_components(), &[0, 1, 0, 0, 0]);
        assert_eq!(m.vertex_components(), vec![0, 0, 0, 0, 1, 1, 1]);
        assert_eq!(m.euler_characteristic(), 3);
        assert_eq!(m.genus(), 0);
        let parts = m.split_components();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].face_to_parent, vec![0, 2, 3, 4]);
        assert_eq!(parts[0].vertex_to_parent, vec![0, 1, 2, 3]);
        assert_eq!(parts[0].mesh.faces(), &t[..]);
        assert_eq!(parts[1].face_to_parent, vec![1]);
        assert_eq!(parts[1].vertex_to_parent, vec![4, 5, 6]);
        assert_eq!(parts[1].mesh.faces(), &[[0, 1, 2]]);

        let e = m.submesh(&[0, 9]).unwrap_err();
        assert_eq!(
            e.to_string(),
            "face selection entry 1 is face 9, but the mesh has 5 faces"
        );
        let e = m.submesh(&[2, 2]).unwrap_err();
        assert_eq!(
            e.to_string(),
            "face selection entry 1 lists face 2 a second time"
        );
        // Two faces of the tetrahedron that share an edge form a valid disk.
        let s = m.submesh(&[2, 0]).unwrap();
        assert_eq!(s.vertex_to_parent, vec![0, 1, 2, 3]);
        assert_eq!(s.mesh.faces(), &[t[1], t[0]]);
        assert_eq!(s.mesh.boundary_loops().len(), 1);
    }
}
