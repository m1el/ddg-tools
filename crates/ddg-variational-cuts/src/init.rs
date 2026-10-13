//! Initial level-set functions and the visibility penalty (SPEC §5).
//!
//! Every function here takes and returns lengths in the unit of the `geometry` it is given, which
//! is the unit [`Flow::new`](crate::Flow::new) expects. A level-set function is a start state like
//! any other: [`Flow::new`](crate::Flow::new) takes any finite values as they are, and
//! [`redistanced`] turns them into distances first.

use ddg_core::distance::upwind_distance;
use ddg_core::{Geometry, TriMesh};

use crate::cut::lone_corner;
use crate::error::CutError;
use crate::flow;
use crate::problem::check_length;

/// A geodesic disk around a vertex: `φ_i = dist(center, i) − radius`, so region 0 is the disk.
///
/// The distance is the fast-marching distance of [`ddg_core::distance::upwind_distance`]: first
/// order, never below the true distance, so the disk is slightly smaller than the geodesic one.
/// A vertex at exactly `radius` has `φ = 0` and belongs to the disk (region 0). A vertex in
/// another connected component than `center` is outside: it gets the largest finite value of
/// `φ`, or `radius` if that is larger. If no vertex other than those of the disk exists, the
/// function has no cut; nothing reports that here (the flow does).
///
/// # Errors
/// [`CutError::Geometry`] if `geometry` belongs to another mesh, [`CutError::IndexOutOfRange`]
/// for `center`, [`CutError::InvalidOptions`] unless `radius` is positive and finite.
pub fn geodesic_circle(
    mesh: &TriMesh,
    geometry: &Geometry,
    center: usize,
    radius: f64,
) -> Result<Vec<f64>, CutError> {
    geometry.check_mesh(mesh)?;
    if center >= mesh.n_vertices() {
        return Err(CutError::IndexOutOfRange {
            what: "center vertex",
            index: center,
            bound: mesh.n_vertices(),
        });
    }
    if !(radius.is_finite() && radius > 0.0) {
        return Err(CutError::InvalidOptions {
            reason: format!("the radius must be positive and finite, got {radius}"),
        });
    }
    let distance = upwind_distance(mesh, geometry, &[(center as u32, 0.0)], None)?;
    let farthest = distance
        .values
        .iter()
        .filter(|d| d.is_finite())
        .fold(0.0f64, |m, d| m.max(*d));
    let outside = (farthest - radius).max(radius);
    Ok(distance
        .values
        .iter()
        .map(|d| if d.is_finite() { d - radius } else { outside })
        .collect())
}

/// A level-set function whose regions are the given vertex labels: `labels[i] = 0` puts vertex
/// `i` into region 0 (`φ_i < 0`), `1` into region 1 (`φ_i > 0`).
///
/// The cut then crosses exactly the edges whose ends have different labels. In a face with one
/// vertex `a` of one label and two of the other, the segment through the two edge midpoints is
/// parallel to the edge opposite `a`, at half the height `h_a` of `a`; all three vertices are
/// at that distance from its line. An end of a crossed edge gets `∓` the smallest such `h_a/2`
/// over its faces that the cut passes through; the crossing of an edge is its midpoint when its
/// two ends get the same distance (always on a mesh of congruent triangles), and otherwise lies
/// nearer to the end with the smaller one. Every other vertex gets its signed fast-marching
/// distance from those ([`redistanced`]). A vertex that no front reaches (no vertex of the
/// other label in its component) gets `∓ √(total area)`. Labels that are all equal give a
/// function without a cut.
///
/// # Errors
/// [`CutError::Geometry`], [`CutError::WrongLength`], and [`CutError::InvalidOptions`] for a
/// label other than 0 or 1.
pub fn from_vertex_labels(
    mesh: &TriMesh,
    geometry: &Geometry,
    labels: &[u32],
) -> Result<Vec<f64>, CutError> {
    geometry.check_mesh(mesh)?;
    check_length("vertex labels", "vertices", labels.len(), mesh.n_vertices())?;
    if let Some(v) = labels.iter().position(|&l| l > 1) {
        return Err(CutError::InvalidOptions {
            reason: format!(
                "vertex {v} has the label {}; two regions are supported (labels 0 and 1)",
                labels[v]
            ),
        });
    }
    let far = geometry.total_area().sqrt();
    let mut distance = vec![f64::INFINITY; labels.len()];
    let lengths = geometry.edge_lengths();
    for (f, face) in mesh.faces().iter().enumerate() {
        let regions = face.map(|v| labels[v as usize]);
        if let Ok(Some(k)) = lone_corner(regions) {
            // The edge opposite corner k is halfedge 3 f + (k + 1) % 3.
            let opposite = lengths[mesh.edge(3 * f + (k + 1) % 3)];
            let half_height = geometry.face_areas()[f] / opposite;
            for &v in face {
                distance[v as usize] = distance[v as usize].min(half_height);
            }
        }
    }
    let phi: Vec<f64> = labels
        .iter()
        .zip(&distance)
        .map(|(&label, &d)| {
            let d = if d.is_finite() { d } else { far };
            if label == 0 { -d } else { d }
        })
        .collect();
    Ok(flow::redistance(mesh, geometry, &phi)?.0)
}

/// Vertex labels from two-valued face labels: a vertex takes the label whose faces cover the
/// larger angle around it; if the two angles agree to a relative `10⁻⁹` (a straight interface
/// through a regular mesh), label 0.
///
/// A vertex all of whose faces have one label gets that label. The vertices on the interface
/// between the two face labels are the ones this rule decides, so the cut of
/// [`from_face_labels`] does not follow the interface edges: it passes through the faces at the
/// interface, within one face of it. A patch of faces can vanish (a single face of label 1
/// among faces of label 0 has no vertex with more than half of its angle in label 1).
///
/// # Errors
/// [`CutError::Geometry`], [`CutError::WrongLength`], and [`CutError::InvalidOptions`] for a
/// label other than 0 or 1.
pub fn vertex_labels_from_faces(
    mesh: &TriMesh,
    geometry: &Geometry,
    face_labels: &[u32],
) -> Result<Vec<u32>, CutError> {
    geometry.check_mesh(mesh)?;
    check_length("face labels", "faces", face_labels.len(), mesh.n_faces())?;
    if let Some(f) = face_labels.iter().position(|&l| l > 1) {
        return Err(CutError::InvalidOptions {
            reason: format!(
                "face {f} has the label {}; two regions are supported (labels 0 and 1)",
                face_labels[f]
            ),
        });
    }
    let mut angle = vec![[0.0f64; 2]; mesh.n_vertices()];
    for (f, face) in mesh.faces().iter().enumerate() {
        for k in 0..3 {
            angle[face[k] as usize][face_labels[f] as usize] += geometry.corner_angles()[f][k];
        }
    }
    Ok(angle
        .iter()
        .map(|a| u32::from(a[1] > a[0] * (1.0 + 1e-9)))
        .collect())
}

/// A level-set function from two-valued face labels: [`vertex_labels_from_faces`], then
/// [`from_vertex_labels`]. See there for how the vertices on the interface are assigned and
/// where the cut runs.
///
/// # Errors
/// As [`vertex_labels_from_faces`].
pub fn from_face_labels(
    mesh: &TriMesh,
    geometry: &Geometry,
    face_labels: &[u32],
) -> Result<Vec<f64>, CutError> {
    let labels = vertex_labels_from_faces(mesh, geometry, face_labels)?;
    from_vertex_labels(mesh, geometry, &labels)
}

/// Makes an arbitrary function a signed distance to its own cut, without moving the cut.
///
/// 1. The function is divided by the mean of `|∇φ|` over its cut: in a crossed face the linear
///    interpolant has one gradient, and the mean is the harmonic one weighted by the length of
///    the cut in the face, `Σ |pq| / Σ (|pq| / |∇φ|)`. A constant factor does not move any
///    crossing.
/// 2. The redistancing of SPEC §3.12: the ends of crossed edges keep their values; every other
///    vertex gets `∓` its fast-marching distance from them within its own region.
///
/// Returns the new function and the number of vertices no front reached (they keep their
/// value, divided like the rest). A function without a cut is returned unchanged.
///
/// After this, `|∇φ|` across the cut is 1 on average and exactly 1 if the input is a multiple
/// of a signed distance; where the slope of the input varies along the cut, that variation
/// stays (removing it would move the cut). The default flow works with the mean slope
/// ([`SpeedScaling::MeanGradient`](crate::SpeedScaling::MeanGradient)); the stages of the paper
/// assume it is 1. To start from a function as it is, pass it to
/// [`Flow::new`](crate::Flow::new) directly.
///
/// # Errors
/// [`CutError::Geometry`], [`CutError::WrongLength`], [`CutError::NonFiniteLevelSet`].
pub fn redistanced(
    mesh: &TriMesh,
    geometry: &Geometry,
    phi: &[f64],
) -> Result<(Vec<f64>, usize), CutError> {
    geometry.check_mesh(mesh)?;
    check_length(
        "the level-set function",
        "vertices",
        phi.len(),
        mesh.n_vertices(),
    )?;
    if let Some(v) = phi.iter().position(|x| !x.is_finite()) {
        return Err(CutError::NonFiniteLevelSet {
            vertex: v as u32,
            value: phi[v],
        });
    }
    // Per crossed face: the gradient of the interpolant from the cotangent formula
    // |∇φ|² A = ½ Σ cot_k (φ_i − φ_j)², and the length of the segment from the two crossing
    // fractions at the lone corner.
    let lengths = geometry.edge_lengths();
    let (mut length, mut time) = (0.0, 0.0);
    for (f, face) in mesh.faces().iter().enumerate() {
        let values = face.map(|v| phi[v as usize]);
        let Ok(Some(k)) = lone_corner(values.map(crate::cut::region_of)) else {
            continue;
        };
        let cot = geometry.corner_cotans()[f];
        let squared: f64 = (0..3)
            .map(|c| cot[c] * (values[(c + 1) % 3] - values[(c + 2) % 3]).powi(2))
            .sum::<f64>()
            / (2.0 * geometry.face_areas()[f]);
        // Halfedge 3 f + c runs from corner c to corner c + 1.
        let to_next = lengths[mesh.edge(3 * f + k)];
        let to_previous = lengths[mesh.edge(3 * f + (k + 2) % 3)];
        let opposite = lengths[mesh.edge(3 * f + (k + 1) % 3)];
        let s_next = values[k] / (values[k] - values[(k + 1) % 3]);
        let s_previous = values[k] / (values[k] - values[(k + 2) % 3]);
        let dot = 0.5 * (to_next * to_next + to_previous * to_previous - opposite * opposite);
        let segment = ((s_next * to_next).powi(2) + (s_previous * to_previous).powi(2)
            - 2.0 * s_next * s_previous * dot)
            .max(0.0)
            .sqrt();
        let gradient = squared.max(0.0).sqrt();
        if segment > 0.0 && gradient > 0.0 {
            length += segment;
            time += segment / gradient;
        }
    }
    let mean = length / time;
    let scaled: Vec<f64> = if mean.is_finite() && mean > 0.0 {
        phi.iter().map(|x| x / mean).collect()
    } else {
        phi.to_vec()
    };
    let (new, unreached, _) = flow::redistance(mesh, geometry, &scaled)?;
    Ok((new, unreached))
}

/// The visibility penalty of the paper (Sec. 7.1): `η_i = ν_i · max(n_i · e_i, 0)`, with `n_i`
/// the unit normal at vertex `i` (the area-weighted mean of the normals of its faces), `e_i`
/// the unit vector from the vertex to `camera`, and `ν_i ≥ 0` the caller's visibility of the
/// vertex (1 where the camera sees it, 0 where something is in the way; this function casts no
/// rays).
///
/// A vertex facing away from the camera gets 0, as do a vertex at the camera point and a vertex
/// whose faces have no area. The result is a number in `[0, ν_i]` without a unit: scale it to
/// the size of the length weight before storing it in
/// [`CutProblem::penalty`](crate::CutProblem::penalty), which is per unit length of the
/// unit-area surface.
///
/// # Errors
/// [`CutError::WrongLength`], [`CutError::InvalidPenalty`] for a visibility that is negative or
/// not finite, [`CutError::InvalidOptions`] for a camera point or a position that is not
/// finite.
pub fn visibility_penalty(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    camera: [f64; 3],
    visibility: &[f64],
) -> Result<Vec<f64>, CutError> {
    let nv = mesh.n_vertices();
    check_length("positions", "vertices", positions.len(), nv)?;
    check_length("visibility", "vertices", visibility.len(), nv)?;
    if let Some(v) = visibility
        .iter()
        .position(|x| !(x.is_finite() && *x >= 0.0))
    {
        return Err(CutError::InvalidPenalty {
            vertex: v as u32,
            value: visibility[v],
        });
    }
    if camera.iter().any(|x| !x.is_finite()) {
        return Err(CutError::InvalidOptions {
            reason: format!("the camera point {camera:?} is not finite"),
        });
    }
    if let Some(v) = positions
        .iter()
        .position(|p| p.iter().any(|x| !x.is_finite()))
    {
        return Err(CutError::InvalidOptions {
            reason: format!(
                "the position {:?} of vertex {v} is not finite",
                positions[v]
            ),
        });
    }
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mut normal = vec![[0.0f64; 3]; nv];
    for face in mesh.faces() {
        let [a, b, c] = face.map(|v| positions[v as usize]);
        let (u, v) = (sub(b, a), sub(c, a));
        // Twice the area times the unit normal.
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        for &i in face {
            for k in 0..3 {
                normal[i as usize][k] += n[k];
            }
        }
    }
    Ok((0..nv)
        .map(|i| {
            let to_camera = sub(camera, positions[i]);
            let scale = (dot(normal[i], normal[i]) * dot(to_camera, to_camera)).sqrt();
            if scale > 0.0 && scale.is_finite() {
                visibility[i] * (dot(normal[i], to_camera) / scale).max(0.0)
            } else {
                0.0
            }
        })
        .collect())
}
