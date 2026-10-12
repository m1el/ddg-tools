//! Helpers shared by the dataset test and the dataset report example: loading the meshes of
//! the optimal-cones dataset and comparing a layout with the UV map stored in them.

#![allow(dead_code)]

use std::path::PathBuf;

use ddg_bff::quality::fit_similarity;
use ddg_bff::{Cone, discontinuous_edges};
use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};

/// The directory of the dataset (it is not part of the repository).
pub fn dataset_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data")
}

/// The names of the dataset meshes, sorted, or `None` if the dataset is absent.
pub fn dataset_names() -> Option<Vec<String>> {
    let dir = dataset_dir();
    if !dir.is_dir() {
        return None;
    }
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .ok()?
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()
                .and_then(|n| n.strip_suffix(".obj"))
                .map(str::to_string)
        })
        .collect();
    names.sort();
    Some(names)
}

/// One mesh of the dataset.
pub struct DatasetMesh {
    pub name: String,
    pub positions: Vec<[f64; 3]>,
    pub mesh: TriMesh,
    pub geometry: Geometry,
    pub cones: Vec<Cone>,
    /// The UV of every corner as stored in the OBJ file (output of the reference BFF).
    pub reference_uvs: Vec<[f64; 2]>,
    /// The edges across which the stored UV map is discontinuous: the reference's cut, as far
    /// as it can be seen. Where both sides of a cut edge have the same UVs to the six digits of
    /// the file (cones of nearly zero angle), the edge is missing here.
    pub reference_cut: Vec<u32>,
}

impl DatasetMesh {
    /// The reference's cut, completed where it could not be recovered from the stored UVs:
    /// if the recovered edges do not cut the mesh into disks through all cones, the planner
    /// adds the missing edges (keeping all recovered ones). Returns the edges and how many
    /// were added.
    pub fn completed_reference_cut(&self) -> (Vec<u32>, usize) {
        let cone_vertices: Vec<u32> = self.cones.iter().map(|c| c.vertex).collect();
        let plan = ddg_bff::plan_cut(&ddg_bff::CutPlanInput {
            mesh: &self.mesh,
            edge_lengths: self.geometry.edge_lengths(),
            cone_vertices: &cone_vertices,
            may_cut: None,
            must_cut: &self.reference_cut,
        })
        .unwrap_or_else(|e| panic!("{}: completing the stored cut: {e}", self.name));
        let added = plan.cut_edges.len() - self.reference_cut.len();
        (plan.cut_edges, added)
    }
}

pub fn load(name: &str) -> DatasetMesh {
    let dir = dataset_dir();
    let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap_or_else(|e| panic!("{e}"));
    let faces = obj
        .triangles()
        .unwrap_or_else(|| panic!("{name}: not triangles"));
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
    let geometry =
        Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{name}: {e}"));
    let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
        .unwrap_or_else(|e| panic!("{e}"));
    let reference_uvs = obj
        .corner_uvs()
        .unwrap_or_else(|| panic!("{name}: corners without UV"));
    let reference_cut = discontinuous_edges(&mesh, &reference_uvs).unwrap();
    DatasetMesh {
        name: name.to_string(),
        positions: obj.positions,
        mesh,
        geometry,
        cones,
        reference_uvs,
        reference_cut,
    }
}

/// How well one layout matches another after the best similarity per chart.
#[derive(Debug, Clone, Copy)]
pub struct LayoutAgreement {
    /// Largest over the charts: root mean square distance / diagonal of the chart's bounding
    /// box in the reference.
    pub rms: f64,
    /// Largest over the charts: maximum distance / diagonal of the bounding box.
    pub max: f64,
    /// Scale of the similarity of the first chart (reference size / own size).
    pub scale: f64,
}

/// Fits, for every chart, the similarity that maps `mine` onto `reference` (both per corner)
/// and returns the worst relative residuals.
pub fn compare_layouts(
    face_chart: &[u32],
    n_charts: usize,
    mine: &[[f64; 2]],
    reference: &[[f64; 2]],
) -> LayoutAgreement {
    let mut corners: Vec<Vec<usize>> = vec![Vec::new(); n_charts];
    for (f, &c) in face_chart.iter().enumerate() {
        corners[c as usize].extend([3 * f, 3 * f + 1, 3 * f + 2]);
    }
    let mut out = LayoutAgreement {
        rms: 0.0,
        max: 0.0,
        scale: f64::NAN,
    };
    for (chart, list) in corners.iter().enumerate() {
        let p: Vec<[f64; 2]> = list.iter().map(|&c| mine[c]).collect();
        let q: Vec<[f64; 2]> = list.iter().map(|&c| reference[c]).collect();
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for point in &q {
            for k in 0..2 {
                lo[k] = lo[k].min(point[k]);
                hi[k] = hi[k].max(point[k]);
            }
        }
        let diagonal = (hi[0] - lo[0]).hypot(hi[1] - lo[1]);
        let fit = fit_similarity(&p, &q).expect("a chart has corners");
        out.rms = out.rms.max(fit.rms / diagonal);
        out.max = out.max.max(fit.max / diagonal);
        if chart == 0 {
            out.scale = fit.scale;
        }
    }
    out
}

/// Quantiles (minimum, 1 %, median, 99 %, maximum) of the ratio of corresponding edge lengths
/// of two per-corner layouts, after dividing by the median ratio.
pub fn length_ratio_quantiles(mine: &[[f64; 2]], reference: &[[f64; 2]]) -> [f64; 5] {
    let length = |uv: &[[f64; 2]], f: usize, k: usize| {
        let (a, b) = (uv[3 * f + k], uv[3 * f + (k + 1) % 3]);
        (b[0] - a[0]).hypot(b[1] - a[1])
    };
    let mut ratios: Vec<f64> = (0..mine.len() / 3)
        .flat_map(|f| (0..3).map(move |k| (f, k)))
        .map(|(f, k)| length(mine, f, k) / length(reference, f, k))
        .filter(|r| r.is_finite())
        .collect();
    ratios.sort_by(|a, b| a.total_cmp(b));
    let n = ratios.len();
    let at = |q: f64| ratios[((n - 1) as f64 * q).round() as usize];
    let median = at(0.5);
    [at(0.0), at(0.01), median, at(0.99), at(1.0)].map(|r| r / median)
}

/// Total length of a set of edges.
pub fn total_length(geometry: &Geometry, edges: &[u32]) -> f64 {
    edges
        .iter()
        .map(|&e| geometry.edge_lengths()[e as usize])
        .sum()
}

/// Number of triangles with negative orientation in a per-corner layout.
pub fn count_flipped(uv: &[[f64; 2]]) -> usize {
    (0..uv.len() / 3)
        .filter(|&f| {
            ddg_bff::quality::signed_double_area(uv[3 * f], uv[3 * f + 1], uv[3 * f + 2]) < 0.0
        })
        .count()
}

/// A copy of the mesh data with every position coordinate and every cone angle moved by a
/// uniformly random amount of up to half a unit in its sixth significant digit: the size of
/// the rounding that writing six digits to a file causes.
pub fn perturbed_like_file_rounding(data: &DatasetMesh, seed: u64) -> (Geometry, Vec<Cone>) {
    let mut rng = ddg_core::generators::SplitMix64::new(seed);
    let mut jitter = |x: f64| {
        if x == 0.0 {
            return x;
        }
        let unit = 10f64.powi(x.abs().log10().floor() as i32 - 5);
        x + unit * rng.uniform(-0.5, 0.5)
    };
    let positions: Vec<[f64; 3]> = data.positions.iter().map(|p| p.map(&mut jitter)).collect();
    let cones: Vec<Cone> = data
        .cones
        .iter()
        .map(|c| Cone {
            vertex: c.vertex,
            angle: jitter(c.angle),
        })
        .collect();
    let geometry = Geometry::from_positions(&data.mesh, &positions)
        .unwrap_or_else(|e| panic!("{}: perturbed geometry: {e}", data.name));
    (geometry, cones)
}

/// Splits every triangle into four by its edge midpoints, without moving any vertex: the
/// surface stays the same polyhedron.
pub fn subdivide(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut positions = positions.to_vec();
    let mut midpoint: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    let mut mid = |a: u32, b: u32, positions: &mut Vec<[f64; 3]>| {
        *midpoint.entry((a.min(b), a.max(b))).or_insert_with(|| {
            let (p, q) = (positions[a as usize], positions[b as usize]);
            positions.push([
                0.5 * (p[0] + q[0]),
                0.5 * (p[1] + q[1]),
                0.5 * (p[2] + q[2]),
            ]);
            positions.len() as u32 - 1
        })
    };
    let mut out = Vec::with_capacity(4 * faces.len());
    for &[a, b, c] in faces {
        let (ab, bc, ca) = (
            mid(a, b, &mut positions),
            mid(b, c, &mut positions),
            mid(c, a, &mut positions),
        );
        out.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
    }
    (positions, out)
}

/// A regular tetrahedron, octahedron or icosahedron (`n_faces` = 4, 8 or 20), subdivided
/// `level` times with flat faces. Its only curvature is at the original vertices.
pub fn flat_polyhedron(n_faces: usize, level: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let (mut positions, mut faces): (Vec<[f64; 3]>, Vec<[u32; 3]>) = match n_faces {
        4 => (
            vec![
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ],
            vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]],
        ),
        8 => (
            vec![
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            vec![
                [0, 2, 4],
                [2, 1, 4],
                [1, 3, 4],
                [3, 0, 4],
                [2, 0, 5],
                [1, 2, 5],
                [3, 1, 5],
                [0, 3, 5],
            ],
        ),
        20 => {
            let g = ddg_core::generators::icosphere(0).unwrap();
            (g.positions, g.faces)
        }
        _ => panic!("no such polyhedron"),
    };
    for _ in 0..level {
        (positions, faces) = subdivide(&positions, &faces);
    }
    (positions, faces)
}

/// Per vertex of the input mesh, the sum of the signed corner angles of a per-corner layout.
pub fn layout_angle_sums(mesh: &TriMesh, uv: &[[f64; 2]]) -> Vec<f64> {
    let mut sums = vec![0.0; mesh.n_vertices()];
    for (f, face) in mesh.faces().iter().enumerate() {
        for k in 0..3 {
            let (a, b, c) = (
                uv[3 * f + k],
                uv[3 * f + (k + 1) % 3],
                uv[3 * f + (k + 2) % 3],
            );
            let (ux, uy, vx, vy) = (b[0] - a[0], b[1] - a[1], c[0] - a[0], c[1] - a[1]);
            sums[face[k] as usize] += (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
        }
    }
    sums
}

/// The largest relative difference between the length of an edge in a per-corner layout and
/// its length on the surface.
pub fn max_length_deviation(mesh: &TriMesh, geometry: &Geometry, uv: &[[f64; 2]]) -> f64 {
    (0..mesh.n_halfedges())
        .map(|h| {
            let (a, b) = (uv[h], uv[mesh.next(h)]);
            let l = geometry.edge_lengths()[mesh.edge(h)];
            ((b[0] - a[0]).hypot(b[1] - a[1]) / l - 1.0).abs()
        })
        .fold(0.0, f64::max)
}

/// The cone angle that the stored UV map realizes at every listed cone, measured on the
/// boundary polygon of the cut map alone (so that flipped or degenerate triangles do not
/// matter): at the copies of a cone vertex the interior angles of the polygon sum to
/// `2π − µ` modulo `2π`. Returned is the representative closest to the listed angle.
pub fn stored_cone_angles(data: &DatasetMesh) -> Vec<f64> {
    let cut = ddg_bff::cut_mesh(&data.mesh, &data.reference_cut).expect("stored cut");
    let mut vertex_uv = vec![[0.0; 2]; cut.mesh.n_vertices()];
    for (f, face) in cut.mesh.faces().iter().enumerate() {
        for k in 0..3 {
            vertex_uv[face[k] as usize] = data.reference_uvs[3 * f + k];
        }
    }
    let mut interior = vec![0.0; data.mesh.n_vertices()];
    for boundary in cut.mesh.boundary_loops() {
        let n = boundary.len();
        for i in 0..n {
            let (a, b, c) = (
                vertex_uv[boundary[(i + n - 1) % n] as usize],
                vertex_uv[boundary[i] as usize],
                vertex_uv[boundary[(i + 1) % n] as usize],
            );
            let (ux, uy, vx, vy) = (b[0] - a[0], b[1] - a[1], c[0] - b[0], c[1] - b[1]);
            let turn = (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
            interior[cut.vertex_to_original[boundary[i] as usize] as usize] +=
                std::f64::consts::PI - turn;
        }
    }
    let tau = std::f64::consts::TAU;
    data.cones
        .iter()
        .map(|c| {
            let difference = (tau - interior[c.vertex as usize] - c.angle).rem_euclid(tau);
            c.angle
                + if difference > tau / 2.0 {
                    difference - tau
                } else {
                    difference
                }
        })
        .collect()
}

/// The connected sum of a surface with its mirror image: face `f` is removed from the surface
/// and from a copy reflected in the plane of that face, and the two are glued along the three
/// edges of the hole. The copy has the same intrinsic geometry with the opposite orientation,
/// so the result is a consistently oriented surface with exact positions; its genus is twice
/// the genus of the input (a torus gives a surface of genus 2), and it has twice the boundary
/// loops. The copy may pass through the original in space, which no intrinsic method notices.
pub fn mirror_double(
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    f: usize,
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let n = positions.len() as u32;
    let hole = faces[f];
    let [p0, p1, p2] = hole.map(|v| positions[v as usize]);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let (e1, e2) = (sub(p1, p0), sub(p2, p0));
    let mut normal = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    normal = normal.map(|x| x / norm);
    // The copy of vertex v: itself on the hole, a new vertex elsewhere.
    let mut copy = vec![0u32; positions.len()];
    let mut out_positions = positions.to_vec();
    for v in 0..n {
        if hole.contains(&v) {
            copy[v as usize] = v;
        } else {
            copy[v as usize] = out_positions.len() as u32;
            let d = sub(positions[v as usize], p0);
            let height = d[0] * normal[0] + d[1] * normal[1] + d[2] * normal[2];
            let p = positions[v as usize];
            out_positions.push([
                p[0] - 2.0 * height * normal[0],
                p[1] - 2.0 * height * normal[1],
                p[2] - 2.0 * height * normal[2],
            ]);
        }
    }
    let mut out_faces = Vec::with_capacity(2 * faces.len() - 2);
    for (g, face) in faces.iter().enumerate() {
        if g != f {
            out_faces.push(*face);
        }
    }
    for (g, face) in faces.iter().enumerate() {
        if g != f {
            let [a, b, c] = face.map(|v| copy[v as usize]);
            out_faces.push([a, c, b]);
        }
    }
    (out_positions, out_faces)
}

/// Removes `k` faces whose vertices are interior and pairwise different, chosen by a seeded
/// shuffle: every removed face becomes a boundary loop of three edges. Returns the remaining
/// faces (vertices keep their ids), or `None` if fewer than `k` such faces exist.
pub fn punch_holes(
    n_vertices: usize,
    faces: &[[u32; 3]],
    k: usize,
    seed: u64,
) -> Option<Vec<[u32; 3]>> {
    let mesh = TriMesh::new(n_vertices, faces.to_vec()).ok()?;
    let mut rng = ddg_core::generators::SplitMix64::new(seed);
    let mut order: Vec<usize> = (0..faces.len()).collect();
    for i in (1..order.len()).rev() {
        order.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
    }
    // A vertex is blocked if it is on the boundary or next to a vertex of a removed face:
    // the holes then do not touch each other or the boundary.
    let mut blocked: Vec<bool> = (0..n_vertices)
        .map(|v| mesh.is_boundary_vertex(v))
        .collect();
    for v in 0..n_vertices {
        if mesh.is_boundary_vertex(v) {
            for w in mesh.vertex_neighbors(v) {
                blocked[w] = true;
            }
        }
    }
    let mut removed = vec![false; faces.len()];
    let mut n_removed = 0;
    for f in order {
        if n_removed == k {
            break;
        }
        if faces[f].iter().any(|&v| blocked[v as usize]) {
            continue;
        }
        removed[f] = true;
        n_removed += 1;
        for &v in &faces[f] {
            blocked[v as usize] = true;
            for w in mesh.vertex_neighbors(v as usize) {
                blocked[w] = true;
                for x in mesh.vertex_neighbors(w) {
                    blocked[x] = true;
                }
            }
        }
    }
    (n_removed == k).then(|| {
        faces
            .iter()
            .zip(&removed)
            .filter(|&(_, &r)| !r)
            .map(|(f, _)| *f)
            .collect()
    })
}
