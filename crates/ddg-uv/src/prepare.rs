//! Stages 1 and 2: validation, the triangle mesh, the seam cut and the pieces.

use std::time::Instant;

use ddg_bff::{BffError, cut_mesh, mollify_edge_lengths};
use ddg_core::geometry::{Geometry, GeometryError, edge_length};
use ddg_core::topology::{MeshError, NONE, TriMesh, compact_vertices, orient_faces};

use crate::error::{UvError, internal};
use crate::input::MeshInput;
use crate::options::{Degenerate, Orientation, UnwrapOptions};
use crate::report::{DegenerateReport, SeamReport, ValidationReport};

/// A connected piece of the surface after the seams were cut, as a triangle mesh of its own.
///
/// Faces keep their corner slots: corner `k` of face `f` is corner `k` of triangle
/// `face_triangle[f]` of the [`Prepared`] mesh.
#[derive(Debug, Clone)]
pub struct Piece {
    /// The connectivity of the piece.
    pub mesh: TriMesh,
    /// Its metric: the edge lengths of the input (mollified if the report says so).
    pub geometry: Geometry,
    /// For every face, the triangle of the prepared mesh.
    pub face_triangle: Vec<u32>,
    /// For every vertex, the vertex of the caller's mesh it is a copy of.
    pub vertex_input: Vec<u32>,
    /// For every edge, whether it is an edge of a polygon (and may therefore be cut) rather
    /// than a diagonal of the triangulation.
    pub is_polygon_edge: Vec<bool>,
    /// Lowest-numbered polygon of the piece.
    pub first_polygon: u32,
    /// Number of polygons of the piece.
    pub n_polygons: usize,
}

impl Piece {
    /// Whether the piece has no boundary.
    pub fn is_closed(&self) -> bool {
        self.mesh.boundary_loops().is_empty()
    }

    /// The genus of the piece.
    pub fn genus(&self) -> usize {
        self.mesh.genus()
    }

    /// The piece cut along some of its edges (edge ids of `self.mesh`), which must leave it
    /// connected. Faces and corners keep their ids.
    pub(crate) fn cut(&self, edges: &[u32]) -> Result<Piece, UvError> {
        let cut = cut_mesh(&self.mesh, edges).map_err(|e| self.flatten_error(&e))?;
        if cut.mesh.n_components() != 1 {
            return Err(internal(format!(
                "a slit split the piece that contains face {} into {} parts",
                self.first_polygon,
                cut.mesh.n_components()
            )));
        }
        let lengths: Vec<f64> = cut
            .edge_to_original
            .iter()
            .map(|&e| self.geometry.edge_lengths()[e as usize])
            .collect();
        let geometry = Geometry::from_edge_lengths(&cut.mesh, lengths)
            .map_err(|e| internal(format!("the geometry of a cut piece was refused: {e}")))?;
        Ok(Piece {
            geometry,
            face_triangle: self.face_triangle.clone(),
            vertex_input: cut
                .vertex_to_original
                .iter()
                .map(|&v| self.vertex_input[v as usize])
                .collect(),
            is_polygon_edge: cut
                .edge_to_original
                .iter()
                .map(|&e| self.is_polygon_edge[e as usize])
                .collect(),
            first_polygon: self.first_polygon,
            n_polygons: self.n_polygons,
            mesh: cut.mesh,
        })
    }

    /// A flattener error on this piece, with the piece named by its lowest polygon.
    pub(crate) fn flatten_error(&self, error: &BffError) -> UvError {
        UvError::Flatten {
            polygon: self.first_polygon,
            message: error.to_string(),
        }
    }
}

/// The result of stages 1 and 2: the validated triangle mesh of the input, and its pieces
/// after the seams were cut.
#[derive(Debug, Clone)]
pub struct Prepared {
    /// Number of vertices, polygons and polygon corners of the input.
    pub(crate) n_input_vertices: usize,
    pub(crate) n_polygons: usize,
    pub(crate) n_corners: usize,
    /// For every triangle and corner slot, the polygon corner (slots 1 and 2 are exchanged
    /// where a polygon was flipped).
    pub(crate) triangle_corners: Vec<[u32; 3]>,
    /// For every triangle, its polygon.
    pub(crate) triangle_polygon: Vec<u32>,
    /// The pieces, in order of their lowest polygon.
    pub pieces: Vec<Piece>,
    /// The seam edges that were cut, with the extensions of one-edge seams, as vertex pairs of
    /// the caller, `[min, max]`.
    pub(crate) seam_cut_edges: Vec<[u32; 2]>,
    /// The polygon edges as sorted vertex pairs of the caller (for the final checks).
    pub(crate) polygon_edges: Vec<[u32; 2]>,
    /// Per vertex of the caller: 1 = pinned, 2 = prescribed, 4 = forbidden, 0 = none.
    pub(crate) steering_role: Vec<u8>,
    /// Per vertex of the caller: the prescribed cone angle (0 where none is prescribed).
    pub(crate) prescribed_angle: Vec<f64>,
    /// Per vertex of the caller: whether a polygon uses it.
    pub(crate) vertex_used: Vec<bool>,
    /// What validation found and did.
    pub validation: ValidationReport,
    /// Warnings raised so far.
    pub(crate) warnings: Vec<String>,
    /// Seconds spent on validation and on cutting the seams.
    pub(crate) seconds_validate: f64,
    pub(crate) seconds_seams: f64,
}

/// Rewrites a `ddg-core` validation error of the triangulated mesh in the caller's indices.
fn mesh_error(e: &MeshError, vertex_input: &[u32], triangle_polygon: &[u32]) -> UvError {
    let v = |x: u32| vertex_input[x as usize];
    let p = |f: u32| triangle_polygon[f as usize];
    let sorted_edge = |edge: [u32; 2]| {
        let (a, b) = (v(edge[0]), v(edge[1]));
        [a.min(b), a.max(b)]
    };
    match e {
        MeshError::NoFaces => UvError::NoPolygons,
        MeshError::RepeatedVertexInFace {
            face,
            vertices: [a, b, c],
        } => UvError::RepeatedVertexInPolygon {
            polygon: p(*face),
            vertex: v(if a == b || a == c { *a } else { *b }),
        },
        MeshError::DuplicateFace {
            first,
            second,
            vertices,
        } => {
            let (p0, p1) = (p(*first), p(*second));
            let mut vertices = vertices.map(v);
            vertices.sort_unstable();
            UvError::DuplicatePolygons {
                first: p0.min(p1),
                second: p0.max(p1),
                vertices,
            }
        }
        MeshError::NonManifoldEdge { edge, faces } => {
            let mut polygons: Vec<u32> = faces.iter().map(|&f| p(f)).collect();
            polygons.sort_unstable();
            polygons.dedup();
            if polygons.len() == faces.len() {
                UvError::NonManifoldEdge {
                    edge: sorted_edge(*edge),
                    polygons,
                }
            } else {
                // Two triangles of one polygon meet the edge: it is a diagonal there.
                UvError::DiagonalOnEdge {
                    edge: sorted_edge(*edge),
                    polygons,
                    n_triangles: faces.len(),
                }
            }
        }
        MeshError::InconsistentOrientation { edge, faces } => {
            let (p0, p1) = (p(faces[0]), p(faces[1]));
            UvError::InconsistentOrientation {
                edge: sorted_edge(*edge),
                polygons: [p0.min(p1), p0.max(p1)],
            }
        }
        MeshError::NonManifoldVertex { vertex, fans } => UvError::NonManifoldVertex {
            vertex: v(*vertex),
            fans: *fans,
        },
        MeshError::NonOrientable { edge, faces, .. } => {
            let (p0, p1) = (p(faces[0]), p(faces[1]));
            UvError::NonOrientable {
                edge: sorted_edge(*edge),
                polygons: [p0.min(p1), p0.max(p1)],
            }
        }
        // Size limits and index ranges were checked on the input; an unreferenced vertex
        // cannot remain after compaction.
        other => internal(format!("the triangle mesh was refused: {other}")),
    }
}

/// Runs stages 1 and 2. See [`crate::unwrap`] for what they do.
///
/// # Errors
/// The input and mesh variants of [`UvError`].
pub fn prepare(input: &MeshInput<'_>, options: &UnwrapOptions) -> Result<Prepared, UvError> {
    let start = Instant::now();
    options.validate()?;
    input.validate()?;
    let mut warnings: Vec<String> = Vec::new();
    let n_input_vertices = input.positions.len();
    let n_polygons = input.n_polygons();
    let n_corners = input.n_corners();
    let offsets = input.polygon_offsets;
    let corner_vertices = input.corner_vertices;

    // A polygon that uses a vertex twice.
    {
        let mut seen_in = vec![NONE; n_input_vertices];
        for p in 0..n_polygons {
            for &v in &corner_vertices[offsets[p] as usize..offsets[p + 1] as usize] {
                if seen_in[v as usize] == p as u32 {
                    return Err(UvError::RepeatedVertexInPolygon {
                        polygon: p as u32,
                        vertex: v,
                    });
                }
                seen_in[v as usize] = p as u32;
            }
        }
    }
    let mut corner_polygon = vec![0u32; n_corners];
    for p in 0..n_polygons {
        corner_polygon[offsets[p] as usize..offsets[p + 1] as usize].fill(p as u32);
    }
    let triangle_polygon: Vec<u32> = input
        .triangles
        .iter()
        .map(|t| corner_polygon[t[0] as usize])
        .collect();

    // The triangle mesh on the vertices that are used.
    let faces: Vec<[u32; 3]> = input
        .triangles
        .iter()
        .map(|t| t.map(|c| corner_vertices[c as usize]))
        .collect();
    let compact = compact_vertices(n_input_vertices, &faces)
        .map_err(|e| internal(format!("compacting the vertices failed: {e}")))?;
    let vertex_input = compact.new_to_old;
    let input_vertex = compact.old_to_new;
    let vertex_used: Vec<bool> = input_vertex.iter().map(|&v| v != NONE).collect();
    let n_vertices = vertex_input.len();
    let n_loose_vertices = n_input_vertices - n_vertices;

    // Orientation.
    let oriented = orient_faces(n_vertices, &compact.faces)
        .map_err(|e| mesh_error(&e, &vertex_input, &triangle_polygon))?;
    let mut n_flipped_polygons = 0usize;
    if oriented.n_flipped > 0 {
        // The triangles of a polygon are wound alike, so they are flipped together.
        let mut polygon_flipped = vec![false; n_polygons];
        for (t, &flipped) in oriented.flipped.iter().enumerate() {
            polygon_flipped[triangle_polygon[t] as usize] |= flipped;
        }
        n_flipped_polygons = polygon_flipped.iter().filter(|&&f| f).count();
        if options.orientation == Orientation::Reject {
            // Let the mesh constructor name the first inconsistent edge.
            return Err(match TriMesh::new(n_vertices, compact.faces) {
                Err(e) => mesh_error(&e, &vertex_input, &triangle_polygon),
                Ok(_) => internal("orient_faces flipped faces of a consistently oriented mesh"),
            });
        }
        let first = polygon_flipped.iter().position(|&f| f).unwrap_or(0);
        warnings.push(format!(
            "the faces were not consistently oriented: {n_flipped_polygons} of {n_polygons} \
             faces were flipped (the first is face {first}); their UVs are mirrored with \
             respect to their own winding"
        ));
    }
    let triangle_corners: Vec<[u32; 3]> = input
        .triangles
        .iter()
        .zip(&oriented.flipped)
        .map(|(&[a, b, c], &flip)| if flip { [a, c, b] } else { [a, b, c] })
        .collect();
    let mesh = TriMesh::new(n_vertices, oriented.faces)
        .map_err(|e| mesh_error(&e, &vertex_input, &triangle_polygon))?;

    // Which edges are polygon edges: a side of a triangle is one iff its two corners are
    // neighbours in their polygon.
    let next_corner = |c: u32| {
        let p = corner_polygon[c as usize] as usize;
        if c + 1 == offsets[p + 1] {
            offsets[p]
        } else {
            c + 1
        }
    };
    let mut is_polygon_edge = vec![false; mesh.n_edges()];
    for (t, corners) in triangle_corners.iter().enumerate() {
        for k in 0..3 {
            let (from, to) = (corners[k], corners[(k + 1) % 3]);
            if next_corner(from) == to || next_corner(to) == from {
                is_polygon_edge[mesh.edge(3 * t + k)] = true;
            }
        }
    }
    let polygon_edges: Vec<[u32; 2]> = (0..mesh.n_edges())
        .filter(|&e| is_polygon_edge[e])
        .map(|e| {
            let [a, b] = mesh.edges()[e];
            let (a, b) = (vertex_input[a as usize], vertex_input[b as usize]);
            [a.min(b), a.max(b)]
        })
        .collect::<Vec<_>>();
    let mut polygon_edges = polygon_edges;
    polygon_edges.sort_unstable();

    // Edge lengths, and what to do about triangles without area.
    let positions = input.positions;
    let mut lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            edge_length(
                positions[vertex_input[a as usize] as usize],
                positions[vertex_input[b as usize] as usize],
            )
        })
        .collect();
    let mut degenerate = DegenerateReport::default();
    let degenerate_error =
        |faces: &[ddg_core::geometry::DegenerateFace], after_mollification: bool| {
            let first = &faces[0];
            UvError::DegenerateTriangles {
                count: faces.len(),
                polygon: triangle_polygon[first.face as usize],
                vertices: first.vertices.map(|v| vertex_input[v as usize]),
                lengths: first.lengths,
                after_mollification,
            }
        };
    let geometry_error = |e: GeometryError, after_mollification: bool| match e {
        GeometryError::DegenerateFaces { faces } => degenerate_error(&faces, after_mollification),
        GeometryError::InvalidEdgeLength { vertices, .. } => UvError::NonFinite {
            what: "positions",
            position: vertex_input[vertices[0] as usize] as usize,
            value: f64::INFINITY,
        },
        other => internal(format!("the geometry was refused: {other}")),
    };
    let mollify = |lengths: &[f64],
                   relative_delta: f64,
                   report: &mut DegenerateReport|
     -> Result<Vec<f64>, UvError> {
        let mollified = mollify_edge_lengths(&mesh, lengths, relative_delta)
            .map_err(|e| internal(format!("mollification failed: {e}")))?;
        report.mollified = mollified.epsilon > 0.0;
        report.epsilon = mollified.epsilon;
        report.delta = mollified.delta;
        report.n_triangles_below_margin = mollified.n_faces_below_margin;
        report.worst_polygon = mollified
            .worst_face
            .map(|f| triangle_polygon[f as usize]);
        Ok(mollified.edge_lengths)
    };
    let geometry = match options.degenerate {
        Degenerate::MollifyAlways { relative_delta } => {
            // Count what the geometry itself would have refused, for the report.
            if let Err(GeometryError::DegenerateFaces { faces }) =
                Geometry::from_edge_lengths(&mesh, lengths.clone())
            {
                degenerate.n_degenerate_triangles = faces.len();
            }
            lengths = mollify(&lengths, relative_delta, &mut degenerate)?;
            Geometry::from_edge_lengths(&mesh, lengths).map_err(|e| geometry_error(e, true))?
        }
        policy => match Geometry::from_edge_lengths(&mesh, lengths.clone()) {
            Ok(geometry) => geometry,
            Err(GeometryError::DegenerateFaces { faces }) => {
                degenerate.n_degenerate_triangles = faces.len();
                match policy {
                    Degenerate::MollifyIfNeeded { relative_delta } => {
                        lengths = mollify(&lengths, relative_delta, &mut degenerate)?;
                        Geometry::from_edge_lengths(&mesh, lengths)
                            .map_err(|e| geometry_error(e, true))?
                    }
                    _ => return Err(degenerate_error(&faces, false)),
                }
            }
            Err(e) => return Err(geometry_error(e, false)),
        },
    };
    if degenerate.mollified {
        warnings.push(format!(
            "{} triangles of the triangulation have no area ({} miss the margin {:.3e}); all \
             edge lengths were lengthened by {:.3e} (intrinsic mollification), the worst \
             triangle is in face {}",
            degenerate.n_degenerate_triangles,
            degenerate.n_triangles_below_margin,
            degenerate.delta,
            degenerate.epsilon,
            degenerate.worst_polygon.unwrap_or(0),
        ));
    }
    let seconds_validate = start.elapsed().as_secs_f64();

    // ----- Stage 2: seams -----
    let start = Instant::now();
    let mut seams = SeamReport {
        n_given: 0,
        ..SeamReport::default()
    };
    let mut is_seam = vec![false; mesh.n_edges()];
    {
        let mut given: Vec<[u32; 2]> = input
            .seam_edges
            .iter()
            .map(|&[a, b]| [a.min(b), a.max(b)])
            .collect();
        given.sort_unstable();
        given.dedup();
        seams.n_given = given.len();
        let mut first_ignored: Option<[u32; 2]> = None;
        for &[a, b] in &given {
            let (ca, cb) = (input_vertex[a as usize], input_vertex[b as usize]);
            let edge = if ca == NONE || cb == NONE {
                None
            } else {
                mesh.find_edge(ca as usize, cb as usize)
            };
            match edge {
                Some(e) if is_polygon_edge[e] => {
                    if mesh.is_boundary_edge(e) {
                        seams.n_on_boundary += 1;
                    } else {
                        is_seam[e] = true;
                        seams.n_cut += 1;
                    }
                }
                _ => {
                    seams.n_not_edges += 1;
                    first_ignored.get_or_insert([a, b]);
                }
            }
        }
        if let Some([a, b]) = first_ignored {
            warnings.push(format!(
                "{} seam edges are not edges of a face and were ignored (the first is \
                 ({a}, {b}))",
                seams.n_not_edges
            ));
        }
    }
    // A seam of one edge between two interior vertices that no other seam touches cannot be
    // represented by an indexed mesh (its two sides would be the same vertex pair). It is
    // extended by the shortest polygon edge at one of its ends.
    let mut seam_edges: Vec<u32> = (0..mesh.n_edges() as u32)
        .filter(|&e| is_seam[e as usize])
        .collect();
    if !seam_edges.is_empty() {
        let mut seam_degree = vec![0u32; n_vertices];
        for &e in &seam_edges {
            for v in mesh.edge_vertices(e as usize) {
                seam_degree[v] += 1;
            }
        }
        let splits = |v: usize, degree: &[u32]| mesh.is_boundary_vertex(v) || degree[v] >= 2;
        let mut first_extended: Option<[u32; 2]> = None;
        for &e in &seam_edges.clone() {
            let [a, b] = mesh.edge_vertices(e as usize);
            if splits(a, &seam_degree) || splits(b, &seam_degree) {
                continue;
            }
            let mut best: Option<(f64, usize)> = None;
            for v in [a, b] {
                for h in mesh.outgoing_halfedges(v) {
                    let candidate = mesh.edge(h);
                    if candidate == e as usize
                        || !is_polygon_edge[candidate]
                        || is_seam[candidate]
                        || mesh.is_boundary_edge(candidate)
                    {
                        continue;
                    }
                    let key = (geometry.edge_lengths()[candidate], candidate);
                    if best.is_none_or(|b| key < b) {
                        best = Some(key);
                    }
                }
            }
            let Some((_, extension)) = best else {
                return Err(internal(format!(
                    "the one-edge seam ({}, {}) cannot be extended",
                    vertex_input[a], vertex_input[b]
                )));
            };
            is_seam[extension] = true;
            for v in mesh.edge_vertices(extension) {
                seam_degree[v] += 1;
            }
            seams.n_extension_edges += 1;
            first_extended.get_or_insert([vertex_input[a], vertex_input[b]]);
        }
        if let Some([a, b]) = first_extended {
            warnings.push(format!(
                "{} seams consist of a single edge between two interior vertices; such a slit \
                 has no two sides, so each was extended by the shortest neighbouring edge (the \
                 first is ({}, {}))",
                seams.n_extension_edges,
                a.min(b),
                a.max(b)
            ));
            seam_edges = (0..mesh.n_edges() as u32)
                .filter(|&e| is_seam[e as usize])
                .collect();
        }
    }
    let mut seam_cut_edges: Vec<[u32; 2]> = seam_edges
        .iter()
        .map(|&e| {
            let [a, b] = mesh.edge_vertices(e as usize);
            let (a, b) = (vertex_input[a], vertex_input[b]);
            [a.min(b), a.max(b)]
        })
        .collect();
    seam_cut_edges.sort_unstable();

    // The pieces.
    let whole = Piece {
        mesh,
        geometry,
        face_triangle: (0..triangle_polygon.len() as u32).collect(),
        vertex_input,
        is_polygon_edge,
        first_polygon: 0,
        n_polygons,
    };
    let mut pieces: Vec<Piece> = if seam_edges.is_empty() && whole.mesh.n_components() == 1 {
        vec![whole]
    } else {
        split_pieces(&whole, &seam_edges, &triangle_polygon, n_polygons)?
    };
    pieces.sort_by_key(|piece| piece.first_polygon);
    seams.n_pieces = pieces.len();

    let mut steering_role = vec![0u8; n_input_vertices];
    let mut prescribed_angle = vec![0.0; n_input_vertices];
    for &v in input.pinned_vertices {
        steering_role[v as usize] |= 1;
    }
    for cone in input.prescribed_cones {
        steering_role[cone.vertex as usize] |= 2;
        prescribed_angle[cone.vertex as usize] = cone.angle;
    }
    for &v in input.forbidden_vertices {
        steering_role[v as usize] |= 4;
    }

    Ok(Prepared {
        n_input_vertices,
        n_polygons,
        n_corners,
        triangle_corners,
        triangle_polygon,
        pieces,
        seam_cut_edges,
        polygon_edges,
        steering_role,
        prescribed_angle,
        vertex_used,
        validation: ValidationReport {
            n_loose_vertices,
            n_flipped_polygons,
            degenerate,
            seams,
        },
        warnings,
        seconds_validate,
        seconds_seams: start.elapsed().as_secs_f64(),
    })
}

/// Cuts `whole` along `seam_edges` and returns the connected components as pieces.
fn split_pieces(
    whole: &Piece,
    seam_edges: &[u32],
    triangle_polygon: &[u32],
    n_polygons: usize,
) -> Result<Vec<Piece>, UvError> {
    let cut = cut_mesh(&whole.mesh, seam_edges)
        .map_err(|e| internal(format!("cutting the seams failed: {e}")))?;
    let lengths = whole.geometry.edge_lengths();
    let mut polygon_seen = vec![false; n_polygons];
    let mut pieces = Vec::with_capacity(cut.mesh.n_components());
    for sub in cut.mesh.split_components() {
        let mesh = sub.mesh;
        // Edge data through the halfedges, which keep their corner slots.
        let mut piece_lengths = vec![0.0; mesh.n_edges()];
        let mut is_polygon_edge = vec![false; mesh.n_edges()];
        for (f, &parent) in sub.face_to_parent.iter().enumerate() {
            for k in 0..3 {
                let e = mesh.edge(3 * f + k);
                let original = whole.mesh.edge(3 * parent as usize + k);
                piece_lengths[e] = lengths[original];
                is_polygon_edge[e] = whole.is_polygon_edge[original];
            }
        }
        let geometry = Geometry::from_edge_lengths(&mesh, piece_lengths)
            .map_err(|e| internal(format!("the geometry of a piece was refused: {e}")))?;
        let face_triangle: Vec<u32> = sub
            .face_to_parent
            .iter()
            .map(|&f| whole.face_triangle[f as usize])
            .collect();
        let mut first_polygon = u32::MAX;
        let mut count = 0usize;
        for &t in &face_triangle {
            let p = triangle_polygon[t as usize];
            first_polygon = first_polygon.min(p);
            if !polygon_seen[p as usize] {
                polygon_seen[p as usize] = true;
                count += 1;
            }
        }
        pieces.push(Piece {
            geometry,
            face_triangle,
            vertex_input: sub
                .vertex_to_parent
                .iter()
                .map(|&v| whole.vertex_input[cut.vertex_to_original[v as usize] as usize])
                .collect(),
            is_polygon_edge,
            first_polygon,
            n_polygons: count,
            mesh,
        });
    }
    Ok(pieces)
}
