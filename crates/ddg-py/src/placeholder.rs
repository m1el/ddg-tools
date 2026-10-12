//! **STAND-IN PIPELINE. This is not the unwrap.**
//!
//! Until `ddg-uv` exists, [`planar_projection`] takes its place behind
//! [`crate::pipeline::run_pipeline`] so that the Python binding and the Blender add-on can be
//! built and tested end to end. It has the pipeline's signature and honours the same contract
//! for its output, but what it computes is deliberately trivial:
//!
//! 1. The triangulated mesh is validated with `ddg_core::TriMesh::new` (after dropping vertices
//!    that no face uses), so a non-manifold or inconsistently oriented mesh is rejected with
//!    ddg-core's message, rewritten in the caller's vertex and polygon indices.
//! 2. Polygons are grouped into charts: two polygons are in one chart when they share an edge
//!    that is not a seam.
//! 3. Every chart is projected orthogonally onto a plane: the plane perpendicular to the chart's
//!    area-weighted normal, or, for a chart whose normals cancel (a closed surface), the
//!    coordinate plane perpendicular to the thinnest side of its bounding box.
//!
//! There is no cone optimization and no conformal flattening. `lambda` and
//! `positive_cones_only` have no effect. The "cones" returned are the steering echoed back
//! (pinned vertices with angle 0, prescribed cones with their angle), and the "cut edges" are the
//! seams that were honoured plus the edges where the projection folds over (a front-facing
//! polygon meets a back-facing one). The report says `placeholder: true` and carries a warning.
//!
//! Delete this file when `ddg-uv` is wired in.

use ddg_core::topology::{MeshError, TriMesh, compact_vertices};

use crate::contract::{Cone, Report, UnwrapInput, UnwrapOptions, UnwrapOutput};
use crate::error::UnwrapError;

/// The warning every placeholder result carries.
pub const PLACEHOLDER_WARNING: &str = "ddg_tools is built with a stand-in pipeline: the UVs are a \
     planar projection per chart, not a conformal map, and no cones were optimized";

type V3 = [f64; 3];

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn edge_key(a: u32, b: u32) -> u64 {
    (u64::from(a.min(b)) << 32) | u64::from(a.max(b))
}
#[inline]
fn key_edge(key: u64) -> [u32; 2] {
    [(key >> 32) as u32, key as u32]
}

/// Rewrites a ddg-core validation error in the caller's indices: `new_to_old` maps the compacted
/// vertex ids back, `triangle_polygon` maps triangles to polygons ("faces" for the user).
fn describe_mesh_error(e: &MeshError, new_to_old: &[u32], triangle_polygon: &[u32]) -> String {
    let v = |x: u32| new_to_old[x as usize];
    let p = |f: u32| triangle_polygon[f as usize];
    match e {
        MeshError::RepeatedVertexInFace {
            face,
            vertices: [a, b, c],
        } => {
            let repeated = if a == b || a == c { a } else { b };
            format!(
                "face {} uses vertex {} more than once",
                p(*face),
                v(*repeated)
            )
        }
        MeshError::DuplicateFace {
            first,
            second,
            vertices: [a, b, c],
        } => {
            let (p0, p1) = (p(*first), p(*second));
            let (a, b, c) = (v(*a), v(*b), v(*c));
            if p0 == p1 {
                format!("face {p0} covers the triangle of vertices ({a}, {b}, {c}) twice")
            } else {
                format!(
                    "faces {} and {} are duplicates or overlap: both cover the triangle of \
                     vertices ({a}, {b}, {c})",
                    p0.min(p1),
                    p0.max(p1)
                )
            }
        }
        MeshError::NonManifoldEdge { edge, faces } => {
            let edge = [v(edge[0]).min(v(edge[1])), v(edge[0]).max(v(edge[1]))];
            let mut polygons: Vec<u32> = faces.iter().map(|&f| p(f)).collect();
            polygons.sort_unstable();
            polygons.dedup();
            if polygons.len() == faces.len() {
                // The usual case: every triangle on the edge comes from a different polygon.
                // Use ddg-core's own wording.
                MeshError::NonManifoldEdge {
                    edge,
                    faces: polygons,
                }
                .to_string()
            } else {
                let list: Vec<String> = polygons.iter().map(|p| p.to_string()).collect();
                format!(
                    "edge ({}, {}) is shared by {} triangles of the triangulation, in faces {}: \
                     a diagonal of a face coincides with another edge",
                    edge[0],
                    edge[1],
                    faces.len(),
                    list.join(", ")
                )
            }
        }
        MeshError::InconsistentOrientation { edge, faces } => {
            let (a, b) = (v(edge[0]).min(v(edge[1])), v(edge[0]).max(v(edge[1])));
            let (p0, p1) = (p(faces[0]).min(p(faces[1])), p(faces[0]).max(p(faces[1])));
            format!(
                "faces {p0} and {p1} traverse their shared edge ({a}, {b}) in the same \
                 direction: the mesh is not consistently oriented (one of them is flipped)"
            )
        }
        MeshError::NonManifoldVertex { vertex, fans } => MeshError::NonManifoldVertex {
            vertex: v(*vertex),
            fans: *fans,
        }
        .to_string(),
        other => other.to_string(),
    }
}

/// Union-find over polygons.
struct Sets(Vec<u32>);

impl Sets {
    fn new(n: usize) -> Sets {
        Sets((0..n as u32).collect())
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.0[x as usize] != x {
            let parent = self.0[x as usize];
            self.0[x as usize] = self.0[parent as usize];
            x = parent;
        }
        x
    }
    /// Joins two sets; the lower root wins, so a set's root is its lowest member.
    fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0[ra.max(rb) as usize] = ra.min(rb);
        }
    }
}

/// The stand-in for the real pipeline: one planar projection per chart.
///
/// See the [module documentation](self). The output satisfies the contract of
/// [`UnwrapOutput`]; it is deterministic (no hashing, no threads).
///
/// # Errors
/// [`UnwrapError::Mesh`] when the triangulated mesh is not an oriented manifold.
pub fn planar_projection(
    input: &UnwrapInput,
    _options: &UnwrapOptions,
) -> Result<UnwrapOutput, UnwrapError> {
    let positions = input.positions();
    let corner_vertices = input.corner_vertices();
    let n_polygons = input.n_polygons();
    let corner_polygon = input.corner_polygons();
    let triangle_polygon: Vec<u32> = input
        .triangles()
        .iter()
        .map(|t| corner_polygon[t[0] as usize])
        .collect();

    // 1. Validate through ddg-core. Vertices that no face uses do not matter for a UV map, so
    //    they are removed first (explicitly, and counted in the report).
    let faces: Vec<[u32; 3]> = input
        .triangles()
        .iter()
        .map(|t| t.map(|c| corner_vertices[c as usize]))
        .collect();
    let compacted = compact_vertices(input.n_vertices(), &faces)
        .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
    let n_unreferenced = input.n_vertices() - compacted.new_to_old.len();
    let new_to_old = compacted.new_to_old;
    TriMesh::new(new_to_old.len(), compacted.faces)
        .map_err(|e| UnwrapError::Mesh(describe_mesh_error(&e, &new_to_old, &triangle_polygon)))?;

    // 2. Charts: polygons joined across the edges that are not seams.
    let mut sides: Vec<(u64, u32)> = Vec::with_capacity(input.n_corners());
    for p in 0..n_polygons {
        let range = input.polygon_corners(p);
        for c in range.clone() {
            let next = if c + 1 == range.end {
                range.start
            } else {
                c + 1
            };
            sides.push((
                edge_key(corner_vertices[c], corner_vertices[next]),
                p as u32,
            ));
        }
    }
    sides.sort_unstable();
    let seams = input.seam_edges();
    let mut seam_used = vec![false; seams.len()];
    let mut sets = Sets::new(n_polygons);
    // Interior edges as (key, polygon, polygon, is seam).
    let mut interior: Vec<(u64, u32, u32, bool)> = Vec::new();
    let mut i = 0;
    while i < sides.len() {
        let key = sides[i].0;
        let mut j = i + 1;
        while j < sides.len() && sides[j].0 == key {
            j += 1;
        }
        let seam = seams.binary_search(&key_edge(key)).ok();
        if let Some(s) = seam {
            seam_used[s] = true;
        }
        match j - i {
            1 => {}
            2 => {
                let (p, q) = (sides[i].1, sides[i + 1].1);
                if seam.is_none() {
                    sets.union(p, q);
                }
                interior.push((key, p, q, seam.is_some()));
            }
            n => {
                // Unreachable after the validation above; reported instead of trusted.
                let [a, b] = key_edge(key);
                return Err(UnwrapError::Mesh(format!(
                    "edge ({a}, {b}) is used by {n} polygon sides"
                )));
            }
        }
        i = j;
    }
    let n_seams_ignored = seam_used.iter().filter(|&&used| !used).count();
    // Number the charts by their lowest polygon.
    let mut polygon_charts = vec![0u32; n_polygons];
    let mut n_charts = 0u32;
    for p in 0..n_polygons {
        let root = sets.find(p as u32) as usize;
        if root == p {
            polygon_charts[p] = n_charts;
            n_charts += 1;
        } else {
            polygon_charts[p] = polygon_charts[root];
        }
    }
    let n_charts = n_charts as usize;

    // 3. One projection plane per chart.
    let mut polygon_normal = vec![[0.0f64; 3]; n_polygons];
    let mut chart_normal = vec![[0.0f64; 3]; n_charts];
    let mut chart_area2 = vec![0.0f64; n_charts];
    for (t, tri) in input.triangles().iter().enumerate() {
        let [a, b, c] = tri.map(|corner| positions[corner_vertices[corner as usize] as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let p = triangle_polygon[t] as usize;
        polygon_normal[p] = add(polygon_normal[p], n);
        let chart = polygon_charts[p] as usize;
        chart_normal[chart] = add(chart_normal[chart], n);
        chart_area2[chart] += norm(n);
    }
    let mut lo = vec![[f64::INFINITY; 3]; n_charts];
    let mut hi = vec![[f64::NEG_INFINITY; 3]; n_charts];
    for (c, &vertex) in corner_vertices.iter().enumerate() {
        let chart = polygon_charts[corner_polygon[c] as usize] as usize;
        let x = positions[vertex as usize];
        for k in 0..3 {
            lo[chart][k] = lo[chart][k].min(x[k]);
            hi[chart][k] = hi[chart][k].max(x[k]);
        }
    }
    // (axis, e1, e2) per chart, right-handed.
    let frames: Vec<[V3; 3]> = (0..n_charts)
        .map(|chart| {
            let n = chart_normal[chart];
            let length = norm(n);
            let axis = if length > 1e-6 * chart_area2[chart] {
                scale(n, 1.0 / length)
            } else {
                // The normals cancel: use the thinnest side of the bounding box; z, y, x on ties.
                let mut best = 2;
                for k in [1, 0] {
                    if hi[chart][k] - lo[chart][k] < hi[chart][best] - lo[chart][best] {
                        best = k;
                    }
                }
                let mut axis = [0.0; 3];
                axis[best] = 1.0;
                axis
            };
            // Helper: the coordinate axis least aligned with `axis` (lowest index on ties).
            let mut k = 0;
            for candidate in [1, 2] {
                if axis[candidate].abs() < axis[k].abs() {
                    k = candidate;
                }
            }
            let mut helper = [0.0; 3];
            helper[k] = 1.0;
            let e1 = sub(helper, scale(axis, dot(helper, axis)));
            let e1 = scale(e1, 1.0 / norm(e1));
            let e2 = cross(axis, e1);
            [axis, e1, e2]
        })
        .collect();

    let uvs: Vec<[f64; 2]> = corner_vertices
        .iter()
        .enumerate()
        .map(|(c, &vertex)| {
            let [_, e1, e2] = frames[polygon_charts[corner_polygon[c] as usize] as usize];
            let x = positions[vertex as usize];
            [dot(x, e1), dot(x, e2)]
        })
        .collect();

    // "Cut" edges: honoured seams, and folds of the projection inside a chart.
    let front = |p: u32| {
        dot(
            polygon_normal[p as usize],
            frames[polygon_charts[p as usize] as usize][0],
        ) >= 0.0
    };
    let mut n_folds = 0usize;
    let mut cut_edges: Vec<[u32; 2]> = Vec::new();
    for &(key, p, q, is_seam) in &interior {
        if is_seam {
            cut_edges.push(key_edge(key));
        } else if front(p) != front(q) {
            n_folds += 1;
            cut_edges.push(key_edge(key));
        }
    }
    // `interior` is in key order, so `cut_edges` is sorted and unique already.

    // "Cones": the steering, echoed.
    let mut cones: Vec<Cone> = input
        .pinned_vertices()
        .iter()
        .map(|&vertex| Cone { vertex, angle: 0.0 })
        .chain(input.prescribed_cones().iter().copied())
        .collect();
    cones.sort_by_key(|c| c.vertex);

    let mut report = Report::new();
    report.set("n_unreferenced_vertices", n_unreferenced);
    report.set("n_seam_edges_ignored", n_seams_ignored);
    report.set("n_fold_edges", n_folds);
    report.warn(PLACEHOLDER_WARNING);
    if n_seams_ignored > 0 {
        report.warn(format!(
            "{n_seams_ignored} seam {} not on any face and {} ignored",
            if n_seams_ignored == 1 {
                "edge is"
            } else {
                "edges are"
            },
            if n_seams_ignored == 1 { "was" } else { "were" },
        ));
    }

    Ok(UnwrapOutput {
        uvs,
        polygon_charts,
        cones,
        cut_edges,
        report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contract::test_meshes::*;
    use crate::contract::{MeshBuffers, ReportValue, check_output};

    fn run(buffers: MeshBuffers) -> Result<UnwrapOutput, UnwrapError> {
        let input = UnwrapInput::from_buffers(buffers).unwrap();
        let out = planar_projection(&input, &UnwrapOptions::default())?;
        check_output(&input, &out).unwrap();
        Ok(out)
    }

    #[test]
    fn flat_grid_projects_to_itself() {
        let out = run(grid(3, 2)).unwrap();
        // One chart, facing +z: (u, v) = (x, y).
        assert!(out.polygon_charts.iter().all(|&c| c == 0));
        let input = UnwrapInput::from_buffers(grid(3, 2)).unwrap();
        for (c, uv) in out.uvs.iter().enumerate() {
            let x = input.positions()[input.corner_vertices()[c] as usize];
            assert_eq!(*uv, [x[0], x[1]]);
        }
        assert!(out.cut_edges.is_empty());
        assert!(out.cones.is_empty());
        assert_eq!(out.report.warnings(), &[PLACEHOLDER_WARNING.to_owned()]);
    }

    #[test]
    fn cube_is_one_chart_folded_along_the_bottom() {
        let out = run(cube()).unwrap();
        assert!(out.polygon_charts.iter().all(|&c| c == 0));
        // Closed: the normals cancel, the box is equally thick everywhere, so the axis is z and
        // the fold separates the bottom face (z = 0) from the rest.
        assert_eq!(out.cut_edges, vec![[0, 1], [0, 2], [1, 3], [2, 3]]);
        assert_eq!(out.report.get("n_fold_edges"), Some(&ReportValue::Int(4)));
    }

    #[test]
    fn seams_split_charts_and_are_reported_as_cuts() {
        // Cut the top face (vertices 4, 5, 6, 7) off the cube.
        let mut b = cube();
        b.seam_edges = vec![4, 5, 5, 7, 7, 6, 6, 4];
        let out = run(b).unwrap();
        assert_eq!(out.polygon_charts, vec![0, 1, 0, 0, 0, 0]);
        // The four seams, and no folds: the open box without lid has a net normal (-z), and its
        // side faces are edge-on.
        assert_eq!(out.cut_edges, vec![[4, 5], [4, 6], [5, 7], [6, 7]]);

        // A seam that does not disconnect anything is still a cut; one that is not an edge of
        // the mesh is ignored with a warning.
        let mut b = cube();
        b.seam_edges = vec![4, 5, 0, 7];
        let out = run(b).unwrap();
        assert!(out.polygon_charts.iter().all(|&c| c == 0));
        assert!(out.cut_edges.contains(&[4, 5]));
        assert_eq!(
            out.report.get("n_seam_edges_ignored"),
            Some(&ReportValue::Int(1))
        );
        assert_eq!(
            out.report.warnings()[1],
            "1 seam edge is not on any face and was ignored"
        );
    }

    #[test]
    fn steering_is_echoed_as_cones() {
        let mut b = cube();
        b.pinned_vertices = vec![6, 2];
        b.prescribed_vertices = vec![4];
        b.prescribed_angles = vec![1.25];
        b.forbidden_vertices = vec![0, 1];
        let out = run(b).unwrap();
        assert_eq!(
            out.cones,
            vec![
                Cone {
                    vertex: 2,
                    angle: 0.0
                },
                Cone {
                    vertex: 4,
                    angle: 1.25
                },
                Cone {
                    vertex: 6,
                    angle: 0.0
                },
            ]
        );
    }

    #[test]
    fn non_manifold_edge_is_reported_in_polygon_indices() {
        let e = run(three_fins()).unwrap_err();
        assert_eq!(
            e,
            UnwrapError::Mesh("edge (0, 1) is shared by 3 faces: 0, 1, 2".to_owned())
        );
    }

    #[test]
    fn errors_use_original_vertex_and_polygon_indices() {
        // Two quads glued along (1, 2), then a third on the same edge. Vertex 0 is loose, so the
        // compacted indices differ from the caller's. The quads are fan-triangulated, so quad p
        // is triangles 2p and 2p + 1.
        let positions = [
            [9.0, 9.0, 9.0], // loose
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let b = buffers(&positions, &[&[1, 2, 3, 4], &[2, 1, 5, 6], &[1, 2, 8, 7]]);
        let e = run(b).unwrap_err();
        assert_eq!(
            e,
            UnwrapError::Mesh("edge (1, 2) is shared by 3 faces: 0, 1, 2".to_owned())
        );

        // Two triangles wound the same way across their shared edge.
        let b = buffers(&positions[1..5], &[&[0, 1, 2], &[0, 1, 3]]);
        assert_eq!(
            run(b).unwrap_err(),
            UnwrapError::Mesh(
                "faces 0 and 1 traverse their shared edge (0, 1) in the same direction: the mesh \
                 is not consistently oriented (one of them is flipped)"
                    .to_owned()
            )
        );

        // A bow tie: two triangles that meet in one vertex only (vertex 2 of the caller, after a
        // loose vertex 0).
        let positions = [
            [9.0, 9.0, 9.0],
            [-1.0, -1.0, 0.0],
            [0.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let b = buffers(&positions, &[&[1, 2, 3], &[2, 4, 5]]);
        assert_eq!(
            run(b).unwrap_err(),
            UnwrapError::Mesh(
                "vertex 2 is non-manifold: its incident faces form 2 separate fans".to_owned()
            )
        );

        // The same face twice.
        let b = buffers(&positions, &[&[1, 2, 3], &[3, 1, 2]]);
        assert_eq!(
            run(b).unwrap_err(),
            UnwrapError::Mesh(
                "faces 0 and 1 are duplicates or overlap: both cover the triangle of vertices \
                 (1, 2, 3)"
                    .to_owned()
            )
        );
    }

    #[test]
    fn loose_vertices_are_dropped_and_counted() {
        let mut b = grid(2, 2);
        b.positions.extend([5.0, 5.0, 5.0]);
        let out = run(b).unwrap();
        assert_eq!(
            out.report.get("n_unreferenced_vertices"),
            Some(&ReportValue::Int(1))
        );
    }

    #[test]
    fn tilted_chart_keeps_lengths() {
        // A quad in the plane x = z, i.e. not aligned with any coordinate plane: the projection
        // onto its own plane is an isometry.
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let positions = [[0.0, 0.0, 0.0], [s, 0.0, s], [s, 1.0, s], [0.0, 1.0, 0.0]];
        let out = run(buffers(&positions, &[&[0, 1, 2, 3]])).unwrap();
        let d = |a: [f64; 2], b: [f64; 2]| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
        assert!((d(out.uvs[0], out.uvs[1]) - 1.0).abs() < 1e-12);
        assert!((d(out.uvs[1], out.uvs[2]) - 1.0).abs() < 1e-12);
        assert!((d(out.uvs[0], out.uvs[2]) - 2f64.sqrt()).abs() < 1e-12);
        // Counter-clockwise in UV, like the polygon seen from its normal side.
        let [a, b, c] = [out.uvs[0], out.uvs[1], out.uvs[2]];
        assert!((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0.0);
    }
}
