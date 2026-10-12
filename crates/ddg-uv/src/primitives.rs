//! Polygon meshes shaped like Blender's primitives, for tests, examples and benchmarks:
//! quads, triangle fans and n-gon caps, wound counter-clockwise seen from outside.
//!
//! Every generator returns a [`PolygonMesh`] whose polygons are convex, split into fans (see
//! [`PolygonMesh::from_polygons`]). This is not Blender's tessellation; for convex polygons
//! both are valid triangulations.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use crate::input::PolygonMesh;

/// How one end of a surface of revolution is closed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum End {
    /// A pole vertex on the axis at this height, joined to the nearest ring by triangles.
    Pole(f64),
    /// One polygon with as many corners as the ring has vertices.
    Cap,
    /// Left open: the ring is a boundary loop.
    Open,
}

/// A surface of revolution around the `z` axis: `rings` lists `(radius, height)` from top to
/// bottom, each ring has `segments` vertices, neighbouring rings are joined by quads, and the
/// two ends are closed as `top` and `bottom` say.
///
/// Returns `None` unless there are at least 3 segments, at least one ring (two if both ends
/// are caps or open), and all radii are positive and finite.
pub fn revolution(
    rings: &[(f64, f64)],
    segments: usize,
    top: End,
    bottom: End,
) -> Option<PolygonMesh> {
    let poles = [top, bottom]
        .iter()
        .filter(|e| matches!(e, End::Pole(_)))
        .count();
    if segments < 3
        || rings.is_empty()
        || (rings.len() < 2 && poles == 0)
        || segments.checked_mul(rings.len())? > 4_000_000
        || rings
            .iter()
            .any(|&(r, z)| !(r.is_finite() && r > 0.0 && z.is_finite()))
    {
        return None;
    }
    let s = segments as u32;
    let mut positions = Vec::with_capacity(segments * rings.len() + 2);
    for &(r, z) in rings {
        for j in 0..segments {
            let phi = 2.0 * PI * j as f64 / segments as f64;
            positions.push([r * phi.cos(), r * phi.sin(), z]);
        }
    }
    let ring = |i: usize, j: u32| (i as u32) * s + j % s;
    let mut polygons: Vec<Vec<u32>> = Vec::new();
    match top {
        End::Pole(z) => {
            let pole = positions.len() as u32;
            positions.push([0.0, 0.0, z]);
            for j in 0..s {
                polygons.push(vec![pole, ring(0, j), ring(0, j + 1)]);
            }
        }
        End::Cap => polygons.push((0..s).map(|j| ring(0, j)).collect()),
        End::Open => {}
    }
    for i in 0..rings.len() - 1 {
        for j in 0..s {
            polygons.push(vec![
                ring(i, j),
                ring(i + 1, j),
                ring(i + 1, j + 1),
                ring(i, j + 1),
            ]);
        }
    }
    let last = rings.len() - 1;
    match bottom {
        End::Pole(z) => {
            let pole = positions.len() as u32;
            positions.push([0.0, 0.0, z]);
            for j in 0..s {
                polygons.push(vec![pole, ring(last, j + 1), ring(last, j)]);
            }
        }
        End::Cap => polygons.push((0..s).rev().map(|j| ring(last, j)).collect()),
        End::Open => {}
    }
    Some(PolygonMesh::from_polygons(positions, &polygons))
}

/// Blender's UV sphere: `segments` meridians, `rings` latitude bands, triangle fans at the
/// poles and quads between. The default there is 32 by 16, radius 1.
pub fn uv_sphere(segments: usize, rings: usize, radius: f64) -> Option<PolygonMesh> {
    if rings < 2 || !(radius.is_finite() && radius > 0.0) {
        return None;
    }
    let profile: Vec<(f64, f64)> = (1..rings)
        .map(|i| {
            let theta = PI * i as f64 / rings as f64;
            (radius * theta.sin(), radius * theta.cos())
        })
        .collect();
    revolution(&profile, segments, End::Pole(radius), End::Pole(-radius))
}

/// A capsule: a cylinder of the given `length` between two hemispheres of `rings` bands each.
pub fn capsule(segments: usize, rings: usize, radius: f64, length: f64) -> Option<PolygonMesh> {
    if rings < 1 || !(radius.is_finite() && radius > 0.0 && length.is_finite() && length > 0.0) {
        return None;
    }
    let half = 0.5 * length;
    let mut profile = Vec::with_capacity(2 * rings);
    for i in 1..=rings {
        let theta = 0.5 * PI * i as f64 / rings as f64;
        profile.push((radius * theta.sin(), half + radius * theta.cos()));
    }
    for i in (1..=rings).rev() {
        let theta = 0.5 * PI * i as f64 / rings as f64;
        profile.push((radius * theta.sin(), -half - radius * theta.cos()));
    }
    revolution(
        &profile,
        segments,
        End::Pole(half + radius),
        End::Pole(-half - radius),
    )
}

/// Blender's cylinder with n-gon caps: two polygons of `segments` corners and `segments`
/// quads in each of the `height_segments` bands. The default there is 32 segments, radius 1,
/// depth 2, one band.
pub fn cylinder(
    segments: usize,
    height_segments: usize,
    radius: f64,
    depth: f64,
) -> Option<PolygonMesh> {
    if height_segments < 1 || !(depth.is_finite() && depth > 0.0) {
        return None;
    }
    let profile: Vec<(f64, f64)> = (0..=height_segments)
        .map(|i| (radius, depth * (0.5 - i as f64 / height_segments as f64)))
        .collect();
    revolution(&profile, segments, End::Cap, End::Cap)
}

/// Blender's cone with an n-gon base: `segments` triangles to the apex and one polygon of
/// `segments` corners. The default there is 32 segments, radius 1, depth 2.
pub fn cone(segments: usize, radius: f64, depth: f64) -> Option<PolygonMesh> {
    if !(depth.is_finite() && depth > 0.0) {
        return None;
    }
    revolution(
        &[(radius, -0.5 * depth)],
        segments,
        End::Pole(0.5 * depth),
        End::Cap,
    )
}

/// Blender's torus in quads: `major` segments around the axis, `minor` around the tube. The
/// default there is 48 by 12 with radii 1 and 0.25.
pub fn torus(
    major: usize,
    minor: usize,
    major_radius: f64,
    minor_radius: f64,
) -> Option<PolygonMesh> {
    if major < 3
        || minor < 3
        || major.checked_mul(minor)? > 4_000_000
        || !(minor_radius.is_finite() && minor_radius > 0.0)
        || !(major_radius.is_finite() && major_radius > minor_radius)
    {
        return None;
    }
    let (m, n) = (major as u32, minor as u32);
    let mut positions = Vec::with_capacity(major * minor);
    for i in 0..major {
        let phi = 2.0 * PI * i as f64 / major as f64;
        for j in 0..minor {
            let theta = 2.0 * PI * j as f64 / minor as f64;
            let r = major_radius + minor_radius * theta.cos();
            positions.push([r * phi.cos(), r * phi.sin(), minor_radius * theta.sin()]);
        }
    }
    let at = |i: u32, j: u32| (i % m) * n + j % n;
    let mut polygons = Vec::with_capacity(major * minor);
    for i in 0..m {
        for j in 0..n {
            polygons.push(vec![at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    Some(PolygonMesh::from_polygons(positions, &polygons))
}

/// The boundary of a set of unit cubes ("voxels") as quads: every face of a cube in `cells`
/// that it does not share with another cube of the set, wound outward. Vertices are the
/// lattice points used, scaled by `cell_size`, in lexicographic order.
///
/// The result is a manifold surface if no two cubes of the set touch along an edge or at a
/// corner only; otherwise the pipeline refuses it. Repeated cells count once.
pub fn voxel_surface(cells: &[[i32; 3]], cell_size: f64) -> PolygonMesh {
    let mut occupied: Vec<[i32; 3]> = cells.to_vec();
    occupied.sort_unstable();
    occupied.dedup();
    let inside = |c: [i32; 3]| occupied.binary_search(&c).is_ok();
    // The four corners of each of the six faces, counter-clockwise seen from outside, as
    // offsets from the cell's lowest corner.
    const FACES: [([i32; 3], [[i32; 3]; 4]); 6] = [
        ([-1, 0, 0], [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0]]),
        ([1, 0, 0], [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]]),
        ([0, -1, 0], [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]),
        ([0, 1, 0], [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]]),
        ([0, 0, -1], [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]),
        ([0, 0, 1], [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]),
    ];
    let mut quads: Vec<[[i32; 3]; 4]> = Vec::new();
    for &c in &occupied {
        for (normal, corners) in FACES {
            if !inside([c[0] + normal[0], c[1] + normal[1], c[2] + normal[2]]) {
                quads.push(corners.map(|o| [c[0] + o[0], c[1] + o[1], c[2] + o[2]]));
            }
        }
    }
    let mut index: BTreeMap<[i32; 3], u32> = BTreeMap::new();
    for quad in &quads {
        for &p in quad {
            index.insert(p, 0);
        }
    }
    let mut positions = Vec::with_capacity(index.len());
    for (i, (p, slot)) in index.iter_mut().enumerate() {
        *slot = i as u32;
        positions.push(p.map(|x| f64::from(x) * cell_size));
    }
    let polygons: Vec<Vec<u32>> = quads
        .iter()
        .map(|quad| quad.iter().map(|p| index[p]).collect())
        .collect();
    PolygonMesh::from_polygons(positions, &polygons)
}

/// Blender's cube: six quads, edge length `size` (2 there).
pub fn cube(size: f64) -> PolygonMesh {
    voxel_surface(&[[0, 0, 0]], size)
}

/// A cube of edge length `size` whose faces are grids of `n × n` quads.
pub fn grid_cube(n: usize, size: f64) -> PolygonMesh {
    let n = n.clamp(1, 256) as i32;
    let mut cells = Vec::with_capacity((n * n * n) as usize);
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                cells.push([x, y, z]);
            }
        }
    }
    voxel_surface(&cells, size / f64::from(n))
}

/// A closed surface of genus `holes` in quads: a slab of `(2 · holes + 1) × 3 × 1` blocks of
/// `n³` voxels each, with every second block of the middle row removed.
pub fn holed_slab(holes: usize, n: usize, size: f64) -> PolygonMesh {
    let n = n.clamp(1, 64) as i32;
    let holes = holes.min(64) as i32;
    let mut cells = Vec::new();
    for bx in 0..2 * holes + 1 {
        for by in 0..3 {
            if by == 1 && bx % 2 == 1 {
                continue;
            }
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        cells.push([bx * n + x, by * n + y, z]);
                    }
                }
            }
        }
    }
    voxel_surface(&cells, size / f64::from(n))
}

/// The same mesh with every vertex moved by a pseudo-random amount of at most `amplitude` in
/// each coordinate (deterministic in `seed`).
pub fn jittered(mesh: &PolygonMesh, amplitude: f64, seed: u64) -> PolygonMesh {
    let mut out = mesh.clone();
    ddg_core::generators::perturb(&mut out.positions, [amplitude; 3], seed);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::TriMesh;

    /// Vertices, edges, polygons and the genus of a closed primitive.
    fn counts(mesh: &PolygonMesh) -> (usize, usize, usize, usize, usize) {
        mesh.input().validate().unwrap();
        let faces: Vec<[u32; 3]> = mesh
            .triangles
            .iter()
            .map(|t| t.map(|c| mesh.corner_vertices[c as usize]))
            .collect();
        let tri = TriMesh::new(mesh.positions.len(), faces).unwrap();
        (
            mesh.positions.len(),
            mesh.edges().len(),
            mesh.n_polygons(),
            tri.genus(),
            tri.boundary_loops().len(),
        )
    }

    #[test]
    fn primitives_have_the_expected_topology() {
        assert_eq!(counts(&uv_sphere(32, 16, 1.0).unwrap()), (482, 992, 512, 0, 0));
        assert_eq!(counts(&cube(2.0)), (8, 12, 6, 0, 0));
        assert_eq!(counts(&grid_cube(4, 2.0)), (98, 192, 96, 0, 0));
        assert_eq!(counts(&cylinder(32, 1, 1.0, 2.0).unwrap()), (64, 96, 34, 0, 0));
        assert_eq!(counts(&cone(32, 1.0, 2.0).unwrap()), (33, 64, 33, 0, 0));
        assert_eq!(counts(&capsule(32, 8, 0.5, 1.0).unwrap()), (514, 1056, 544, 0, 0));
        assert_eq!(counts(&torus(48, 12, 1.0, 0.25).unwrap()), (576, 1152, 576, 1, 0));
        let slab = holed_slab(2, 1, 1.0);
        assert_eq!(counts(&slab).3, 2);
        let open = revolution(&[(1.0, 1.0), (1.0, 0.0)], 8, End::Open, End::Cap).unwrap();
        assert_eq!(counts(&open).4, 1);
    }

    #[test]
    fn primitives_are_wound_outward() {
        // The signed volume enclosed by the triangulation is positive.
        for mesh in [
            uv_sphere(12, 6, 1.0).unwrap(),
            cube(2.0),
            cylinder(12, 2, 1.0, 2.0).unwrap(),
            cone(12, 1.0, 2.0).unwrap(),
            capsule(12, 3, 0.5, 1.0).unwrap(),
            torus(12, 6, 1.0, 0.25).unwrap(),
            holed_slab(1, 2, 1.0),
        ] {
            let mut volume = 0.0;
            for t in &mesh.triangles {
                let [a, b, c] = t.map(|k| mesh.positions[mesh.corner_vertices[k as usize] as usize]);
                volume += a[0] * (b[1] * c[2] - b[2] * c[1])
                    + a[1] * (b[2] * c[0] - b[0] * c[2])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]);
            }
            assert!(volume > 0.0, "volume {volume}");
        }
    }

    #[test]
    fn bad_parameters_give_none() {
        assert!(uv_sphere(2, 16, 1.0).is_none());
        assert!(uv_sphere(8, 1, 1.0).is_none());
        assert!(torus(8, 8, 0.2, 0.25).is_none());
        assert!(cone(8, -1.0, 2.0).is_none());
        assert!(cylinder(8, 0, 1.0, 2.0).is_none());
        assert!(capsule(8, 0, 1.0, 2.0).is_none());
    }
}
