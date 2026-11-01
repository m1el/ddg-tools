//! The validated triangle mesh `Surface` (`SPEC.md` §2.1, D1).
//!
//! `ddg_core::TriMesh` requires a consistently oriented mesh; the tangent-point energy, the
//! metric, the barycenter, pins, areas, remeshing and CCD do not depend on the orientation
//! (the kernel uses `|N·v|`), so this crate has its own validated type that accepts
//! non-orientable surfaces (Klein bottle, Möbius band, paper Fig. 26) and reports the
//! orientation class of every connected component.
//!
//! # Validation (in this order)
//!
//! 1. at least one face ([`RsError::NoFaces`]), ids fit in `u32` ([`RsError::TooLarge`]);
//! 2. every corner in range ([`RsError::VertexOutOfRange`]), three distinct vertices per face
//!    ([`RsError::RepeatedVertex`]);
//! 3. no two faces with the same vertex set ([`RsError::DuplicateFace`]);
//! 4. every edge has one or two faces ([`RsError::NonManifoldEdge`], the lowest such edge);
//! 5. the faces around every vertex form one fan, a disk or a half-disk
//!    ([`RsError::NonManifoldVertex`], the lowest such vertex);
//! 6. every vertex belongs to a face ([`RsError::UnreferencedVertex`]).
//!
//! # Conventions
//!
//! * Edges are numbered in lexicographic order of `[min, max]`; edge slot `k` of face `f` is
//!   the edge `(faces[f][k], faces[f][(k + 1) % 3])`.
//! * Components are numbered by their lowest face; [`Component::faces`] and
//!   [`Component::vertices`] are ascending.
//! * The fan of a vertex lists its faces in rotational order. A boundary vertex's fan starts at
//!   the lowest face with a boundary edge at the vertex and ends at the other one; an interior
//!   vertex's fan starts at its lowest face.
//! * A boundary loop starts at its lowest vertex and follows the boundary edge that the edge's
//!   face traverses out of that vertex (on a consistently oriented component, every loop runs
//!   with the faces' orientation, interior to the left); on non-orientable components where
//!   both or neither edge leave the vertex, it takes the edge to the lower neighbour.

use crate::error::RsError;

/// Marker for "no element" in `u32` index arrays (e.g. the missing second face of a boundary
/// edge).
pub const NONE: u32 = u32::MAX;

/// The orientation class of a connected component (`SPEC.md` §2.1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Orientation {
    /// Every interior edge is traversed in opposite directions by its two faces.
    Consistent,
    /// Orientable, but the caller's face order is not consistent: flipping the listed faces
    /// (ascending) makes it consistent. The component's lowest face keeps its orientation. The
    /// caller's order is kept; nothing is flipped.
    Orientable {
        /// The faces to flip, ascending.
        flips: Vec<u32>,
    },
    /// No consistent orientation exists (Klein bottle, Möbius band).
    NonOrientable,
}

/// One connected component of a [`Surface`].
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    faces: Vec<u32>,
    vertices: Vec<u32>,
    n_edges: usize,
    orientation: Orientation,
    inconsistent_edge: Option<[u32; 2]>,
    boundary_loops: Vec<u32>,
}

impl Component {
    /// The faces of the component, ascending.
    pub fn faces(&self) -> &[u32] {
        &self.faces
    }

    /// The vertices of the component, ascending.
    pub fn vertices(&self) -> &[u32] {
        &self.vertices
    }

    /// Number of edges of the component.
    pub fn n_edges(&self) -> usize {
        self.n_edges
    }

    /// The orientation class.
    pub fn orientation(&self) -> &Orientation {
        &self.orientation
    }

    /// The lowest edge (in edge order, as `[min, max]`) whose two faces traverse it in the same
    /// direction in the caller's face order; `None` exactly when the component is
    /// [`Orientation::Consistent`]. Used to locate errors of operations that need an
    /// orientation (signed volume).
    pub fn inconsistent_edge(&self) -> Option<[u32; 2]> {
        self.inconsistent_edge
    }

    /// Indices into [`Surface::boundary_loops`] of the loops of this component, ascending.
    pub fn boundary_loops(&self) -> &[u32] {
        &self.boundary_loops
    }

    /// Euler characteristic `V − E + F` of the component.
    pub fn euler_characteristic(&self) -> i64 {
        self.vertices.len() as i64 - self.n_edges as i64 + self.faces.len() as i64
    }
}

/// A validated manifold triangle mesh; orientation optional (`SPEC.md` §2.1).
#[derive(Debug, Clone, PartialEq)]
pub struct Surface {
    n_vertices: usize,
    faces: Vec<[u32; 3]>,
    edges: Vec<[u32; 2]>,
    face_edges: Vec<[u32; 3]>,
    edge_faces: Vec<[u32; 2]>,
    fan_start: Vec<u32>,
    fan_faces: Vec<u32>,
    vertex_boundary: Vec<bool>,
    face_component: Vec<u32>,
    vertex_component: Vec<u32>,
    components: Vec<Component>,
    boundary_loops: Vec<Vec<u32>>,
}

impl Surface {
    /// Validates the faces (see the module documentation for the checks and their order) and
    /// builds the edges, the edge–face incidence, the fans, the components with their
    /// orientation classes, and the boundary loops.
    ///
    /// # Errors
    /// The topology variants of [`RsError`] listed in the module documentation, each naming the
    /// offending element.
    pub fn new(n_vertices: usize, faces: Vec<[u32; 3]>) -> Result<Surface, RsError> {
        let nf = faces.len();
        if nf == 0 {
            return Err(RsError::NoFaces);
        }
        if n_vertices >= NONE as usize || nf >= (NONE as usize) / 3 {
            return Err(RsError::TooLarge {
                n_vertices,
                n_faces: nf,
            });
        }
        // 2. corners.
        for (f, face) in faces.iter().enumerate() {
            for (k, &v) in face.iter().enumerate() {
                if v as usize >= n_vertices {
                    return Err(RsError::VertexOutOfRange {
                        face: f as u32,
                        corner: k as u8,
                        vertex: v,
                        n_vertices,
                    });
                }
            }
            if face[0] == face[1] || face[1] == face[2] || face[0] == face[2] {
                return Err(RsError::RepeatedVertex {
                    face: f as u32,
                    vertices: *face,
                });
            }
        }
        // 3. duplicates: the pair whose second face is lowest.
        {
            let mut keys: Vec<([u32; 3], u32)> = faces
                .iter()
                .enumerate()
                .map(|(f, t)| {
                    let mut s = *t;
                    s.sort_unstable();
                    (s, f as u32)
                })
                .collect();
            keys.sort_unstable();
            let mut best: Option<(u32, u32, [u32; 3])> = None;
            let mut i = 0;
            while i < keys.len() {
                let mut j = i + 1;
                while j < keys.len() && keys[j].0 == keys[i].0 {
                    j += 1;
                }
                if j - i > 1 {
                    let cand = (keys[i].1, keys[i + 1].1, keys[i].0);
                    if best.is_none_or(|b| cand.1 < b.1) {
                        best = Some(cand);
                    }
                }
                i = j;
            }
            if let Some((first, second, vertices)) = best {
                return Err(RsError::DuplicateFace {
                    first,
                    second,
                    vertices,
                });
            }
        }
        // 4. edges.
        let mut half: Vec<([u32; 2], u32, u8)> = Vec::with_capacity(3 * nf);
        for (f, t) in faces.iter().enumerate() {
            for k in 0..3 {
                let a = t[k];
                let b = t[(k + 1) % 3];
                half.push(([a.min(b), a.max(b)], f as u32, k as u8));
            }
        }
        half.sort_unstable();
        let mut edges: Vec<[u32; 2]> = Vec::new();
        let mut edge_faces: Vec<[u32; 2]> = Vec::new();
        let mut face_edges = vec![[NONE; 3]; nf];
        let mut i = 0;
        while i < half.len() {
            let mut j = i + 1;
            while j < half.len() && half[j].0 == half[i].0 {
                j += 1;
            }
            if j - i > 2 {
                return Err(RsError::NonManifoldEdge {
                    edge: half[i].0,
                    faces: half[i..j].iter().map(|h| h.1).collect(),
                });
            }
            let e = edges.len() as u32;
            edges.push(half[i].0);
            let second = if j - i == 2 { half[i + 1].1 } else { NONE };
            edge_faces.push([half[i].1, second]);
            for h in &half[i..j] {
                face_edges[h.1 as usize][h.2 as usize] = e;
            }
            i = j;
        }
        drop(half);

        // 5. fans.
        let mut degree = vec![0u32; n_vertices + 1];
        for t in &faces {
            for &v in t {
                degree[v as usize + 1] += 1;
            }
        }
        for v in 0..n_vertices {
            degree[v + 1] += degree[v];
        }
        let incident_start = degree;
        let mut fill = incident_start.clone();
        let mut incident = vec![0u32; 3 * nf];
        for (f, t) in faces.iter().enumerate() {
            for &v in t {
                incident[fill[v as usize] as usize] = f as u32;
                fill[v as usize] += 1;
            }
        }
        // 6 (checked after 5, but needs the degrees): unreferenced vertices.
        let mut fan_faces = Vec::with_capacity(3 * nf);
        let mut fan_start = Vec::with_capacity(n_vertices + 1);
        let mut vertex_boundary = vec![false; n_vertices];
        let mut nonmanifold: Option<(u32, u32)> = None;
        let mut unreferenced: Option<(u32, usize)> = None;
        // Neighbours of face `f` around vertex `v`: the other faces of its two edges at `v`.
        let around = |f: u32, v: u32| -> [u32; 2] {
            let t = faces[f as usize];
            let mut out = [NONE; 2];
            let mut n = 0;
            for k in 0..3 {
                if t[k] == v || t[(k + 1) % 3] == v {
                    let e = face_edges[f as usize][k] as usize;
                    let [g0, g1] = edge_faces[e];
                    out[n] = if g0 == f { g1 } else { g0 };
                    n += 1;
                }
            }
            out
        };
        for v in 0..n_vertices {
            fan_start.push(fan_faces.len() as u32);
            let inc = &incident[incident_start[v] as usize..incident_start[v + 1] as usize];
            if inc.is_empty() {
                unreferenced = Some(match unreferenced {
                    None => (v as u32, 1),
                    Some((first, c)) => (first, c + 1),
                });
                continue;
            }
            let mut visited = vec![false; inc.len()];
            let pos = |f: u32| inc.iter().position(|&g| g == f);
            let neighbours: Vec<[u32; 2]> = inc.iter().map(|&f| around(f, v as u32)).collect();
            let mut fans = 0u32;
            let mut first_fan: Vec<u32> = Vec::new();
            let mut first_open = false;
            // Open chains first (from a face with a boundary edge at v), then closed cycles.
            for pass in 0..2 {
                for start in 0..inc.len() {
                    if visited[start] {
                        continue;
                    }
                    let open = neighbours[start].contains(&NONE);
                    if pass == 0 && !open {
                        continue;
                    }
                    fans += 1;
                    let mut chain = Vec::new();
                    let mut prev = NONE;
                    let mut cur = start;
                    loop {
                        visited[cur] = true;
                        chain.push(inc[cur]);
                        let [n0, n1] = neighbours[cur];
                        // For an open start, step away from the boundary side.
                        let next = if prev == NONE {
                            if n0 == NONE { n1 } else { n0 }
                        } else if n0 == prev {
                            n1
                        } else {
                            n0
                        };
                        if next == NONE {
                            break;
                        }
                        let Some(p) = pos(next) else { break };
                        if visited[p] {
                            break;
                        }
                        prev = inc[cur];
                        cur = p;
                    }
                    if fans == 1 {
                        first_fan = chain;
                        first_open = open;
                    }
                }
            }
            if fans > 1 {
                if nonmanifold.is_none() {
                    nonmanifold = Some((v as u32, fans));
                }
                continue;
            }
            vertex_boundary[v] = first_open;
            fan_faces.extend_from_slice(&first_fan);
        }
        fan_start.push(fan_faces.len() as u32);
        if let Some((vertex, fans)) = nonmanifold {
            return Err(RsError::NonManifoldVertex { vertex, fans });
        }
        if let Some((vertex, count)) = unreferenced {
            return Err(RsError::UnreferencedVertex { vertex, count });
        }

        // Components (face BFS in increasing face order: numbered by their lowest face).
        let mut face_component = vec![NONE; nf];
        let mut flip = vec![false; nf];
        let mut comp_faces: Vec<Vec<u32>> = Vec::new();
        let mut conflict: Vec<bool> = Vec::new();
        // Does face f traverse edge e = [min, max] from min to max?
        let forward = |f: u32, e: usize| -> bool {
            let t = faces[f as usize];
            let [a, b] = edges[e];
            (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b)
        };
        for seed in 0..nf {
            if face_component[seed] != NONE {
                continue;
            }
            let c = comp_faces.len() as u32;
            let mut list = vec![seed as u32];
            face_component[seed] = c;
            let mut has_conflict = false;
            let mut head = 0;
            while head < list.len() {
                let f = list[head];
                head += 1;
                for &e in &face_edges[f as usize] {
                    let [g0, g1] = edge_faces[e as usize];
                    let g = if g0 == f { g1 } else { g0 };
                    if g == NONE {
                        continue;
                    }
                    // Consistent: (forward(f) ^ flip f) != (forward(g) ^ flip g).
                    let want =
                        forward(f, e as usize) ^ flip[f as usize] ^ forward(g, e as usize) ^ true;
                    if face_component[g as usize] == NONE {
                        face_component[g as usize] = c;
                        flip[g as usize] = want;
                        list.push(g);
                    } else if flip[g as usize] != want {
                        has_conflict = true;
                    }
                }
            }
            list.sort_unstable();
            comp_faces.push(list);
            conflict.push(has_conflict);
        }
        let n_comp = comp_faces.len();
        let mut vertex_component = vec![NONE; n_vertices];
        for (f, t) in faces.iter().enumerate() {
            for &v in t {
                vertex_component[v as usize] = face_component[f];
            }
        }
        let mut comp_vertices: Vec<Vec<u32>> = vec![Vec::new(); n_comp];
        for v in 0..n_vertices {
            comp_vertices[vertex_component[v] as usize].push(v as u32);
        }
        let mut comp_edges = vec![0usize; n_comp];
        let mut witness: Vec<Option<[u32; 2]>> = vec![None; n_comp];
        for (e, &[f, g]) in edge_faces.iter().enumerate() {
            let c = face_component[f as usize] as usize;
            comp_edges[c] += 1;
            if g != NONE && witness[c].is_none() && forward(f, e) == forward(g, e) {
                witness[c] = Some(edges[e]);
            }
        }

        // Boundary loops: for every boundary vertex, its two boundary edges.
        let mut bnext: Vec<[(u32, u32); 2]> = vec![[(NONE, NONE); 2]; n_vertices];
        for (e, &[_, g]) in edge_faces.iter().enumerate() {
            if g != NONE {
                continue;
            }
            let [a, b] = edges[e];
            for (x, y) in [(a, b), (b, a)] {
                let slot = &mut bnext[x as usize];
                if slot[0].0 == NONE {
                    slot[0] = (y, e as u32);
                } else {
                    slot[1] = (y, e as u32);
                }
            }
        }
        // Does the face of boundary edge e traverse it out of x?
        let leaves = |x: u32, e: u32| -> bool {
            let e = e as usize;
            forward(edge_faces[e][0], e) == (x == edges[e][0])
        };
        let mut boundary_loops: Vec<Vec<u32>> = Vec::new();
        let mut on_loop = vec![false; n_vertices];
        let mut comp_loops: Vec<Vec<u32>> = vec![Vec::new(); n_comp];
        for v in 0..n_vertices as u32 {
            if !vertex_boundary[v as usize] || on_loop[v as usize] {
                continue;
            }
            let [(w0, e0), (w1, e1)] = bnext[v as usize];
            let (l0, l1) = (leaves(v, e0), leaves(v, e1));
            let first = if l0 != l1 {
                if l0 { w0 } else { w1 }
            } else {
                w0.min(w1)
            };
            let mut lp = vec![v];
            on_loop[v as usize] = true;
            let mut prev = v;
            let mut cur = first;
            while cur != v {
                lp.push(cur);
                on_loop[cur as usize] = true;
                let [(n0, _), (n1, _)] = bnext[cur as usize];
                let next = if n0 == prev { n1 } else { n0 };
                prev = cur;
                cur = next;
            }
            comp_loops[vertex_component[v as usize] as usize].push(boundary_loops.len() as u32);
            boundary_loops.push(lp);
        }

        let components = (0..n_comp)
            .map(|c| {
                let orientation = if witness[c].is_none() {
                    Orientation::Consistent
                } else if conflict[c] {
                    Orientation::NonOrientable
                } else {
                    Orientation::Orientable {
                        flips: comp_faces[c]
                            .iter()
                            .copied()
                            .filter(|&f| flip[f as usize])
                            .collect(),
                    }
                };
                Component {
                    faces: std::mem::take(&mut comp_faces[c]),
                    vertices: std::mem::take(&mut comp_vertices[c]),
                    n_edges: comp_edges[c],
                    orientation,
                    inconsistent_edge: witness[c],
                    boundary_loops: std::mem::take(&mut comp_loops[c]),
                }
            })
            .collect();

        Ok(Surface {
            n_vertices,
            faces,
            edges,
            face_edges,
            edge_faces,
            fan_start,
            fan_faces,
            vertex_boundary,
            face_component,
            vertex_component,
            components,
            boundary_loops,
        })
    }

    /// Checks positions against the surface: one per vertex, finite, and every face of
    /// positive finite area (`SPEC.md` §2.3). Nothing is repaired.
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`], [`RsError::NonFinitePosition`] (lowest vertex),
    /// [`RsError::DegenerateFace`] (lowest face).
    pub fn validate_positions(&self, positions: &[[f64; 3]]) -> Result<(), RsError> {
        if positions.len() != self.n_vertices {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        if let Some(v) = positions
            .iter()
            .position(|p| !p.iter().all(|x| x.is_finite()))
        {
            return Err(RsError::NonFinitePosition { vertex: v as u32 });
        }
        for (f, t) in self.faces.iter().enumerate() {
            let p = [
                positions[t[0] as usize],
                positions[t[1] as usize],
                positions[t[2] as usize],
            ];
            let (area, normal, _) = crate::geometry::face_quantities(p);
            if !(area > 0.0 && area.is_finite() && normal.iter().all(|x| x.is_finite())) {
                return Err(RsError::DegenerateFace {
                    face: f as u32,
                    area,
                });
            }
        }
        Ok(())
    }

    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of edges.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// The faces, in the caller's vertex order.
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// The edges as `[min, max]`, in lexicographic order.
    pub fn edges(&self) -> &[[u32; 2]] {
        &self.edges
    }

    /// For every face, its three edge ids; slot `k` is the edge
    /// `(faces[f][k], faces[f][(k + 1) % 3])`.
    pub fn face_edges(&self) -> &[[u32; 3]] {
        &self.face_edges
    }

    /// For every edge, its faces `[lower, higher]`; the second is [`NONE`] on the boundary.
    pub fn edge_faces(&self) -> &[[u32; 2]] {
        &self.edge_faces
    }

    /// Is edge `e` a boundary edge (one face)?
    pub fn is_boundary_edge(&self, e: u32) -> bool {
        self.edge_faces[e as usize][1] == NONE
    }

    /// The id of the edge joining `a` and `b` (either order), if there is one.
    pub fn find_edge(&self, a: u32, b: u32) -> Option<u32> {
        self.edges
            .binary_search(&[a.min(b), a.max(b)])
            .ok()
            .map(|e| e as u32)
    }

    /// The faces around vertex `v` in rotational order (one fan). For a boundary vertex the
    /// first and last faces each have a boundary edge at `v`.
    pub fn fan(&self, v: u32) -> &[u32] {
        let v = v as usize;
        &self.fan_faces[self.fan_start[v] as usize..self.fan_start[v + 1] as usize]
    }

    /// Is `v` on the boundary (its fan is a half-disk)?
    pub fn is_boundary_vertex(&self, v: u32) -> bool {
        self.vertex_boundary[v as usize]
    }

    /// The connected components, numbered by their lowest face.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Number of connected components.
    pub fn n_components(&self) -> usize {
        self.components.len()
    }

    /// The orientation class of component `c` (`SPEC.md` §2.1; returned by reference since
    /// [`Orientation::Orientable`] carries the flip set). Panics if `c >= n_components()`.
    pub fn orientation(&self, c: usize) -> &Orientation {
        &self.components[c].orientation
    }

    /// The component of every face.
    pub fn face_components(&self) -> &[u32] {
        &self.face_component
    }

    /// The component of every vertex.
    pub fn vertex_components(&self) -> &[u32] {
        &self.vertex_component
    }

    /// The boundary loops as vertex cycles (see the module documentation for start and
    /// direction), ordered by their lowest vertex.
    pub fn boundary_loops(&self) -> &[Vec<u32>] {
        &self.boundary_loops
    }

    /// Euler characteristic `V − E + F` of the whole surface.
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices as i64 - self.edges.len() as i64 + self.faces.len() as i64
    }

    /// Is every component [`Orientation::Consistent`]?
    pub fn is_consistently_oriented(&self) -> bool {
        self.components
            .iter()
            .all(|c| c.orientation == Orientation::Consistent)
    }
}
