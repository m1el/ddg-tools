//! Small procedural triangle meshes for tests and benchmarks, one per topology class:
//! sphere ([`icosphere`]), torus ([`torus`]), disk ([`grid_disk`], [`bumpy_grid_disk`]) and
//! annulus ([`annulus`], [`open_cylinder`]), plus a seeded perturbation ([`perturb`]).
//!
//! All generators are deterministic, produce consistently oriented faces (counter-clockwise
//! seen from outside, or from `+z` for the planar ones) and have no dependencies.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt;

use crate::topology::{MeshError, TriMesh};

/// Positions and triangles produced by a generator.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedMesh {
    /// Vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Triangles as vertex indices.
    pub faces: Vec<[u32; 3]>,
}

impl GeneratedMesh {
    /// Builds the validated connectivity. Infallible for unmodified generator output.
    pub fn topology(&self) -> Result<TriMesh, MeshError> {
        TriMesh::new(self.positions.len(), self.faces.clone())
    }
}

/// A generator was called with parameters that do not describe a mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneratorError {
    message: String,
}

impl fmt::Display for GeneratorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for GeneratorError {}

fn bad(message: String) -> GeneratorError {
    GeneratorError { message }
}

fn check_positive(name: &str, value: f64) -> Result<(), GeneratorError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(bad(format!("{name} must be positive and finite, got {value}")))
    }
}

/// The SplitMix64 pseudo-random generator: tiny, fast, and reproducible across platforms.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Creates a generator from a seed. Equal seeds give equal sequences.
    pub fn new(seed: u64) -> Self {
        SplitMix64 { state: seed }
    }

    /// The next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform sample in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// A uniform sample in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// Adds independent uniform noise in `[-amplitude[k], amplitude[k]]` to coordinate `k` of every
/// position. Deterministic in `seed`. Use a zero amplitude to keep an axis fixed (for example
/// `[a, a, 0.0]` keeps a planar mesh planar).
pub fn perturb(positions: &mut [[f64; 3]], amplitude: [f64; 3], seed: u64) {
    let mut rng = SplitMix64::new(seed);
    for p in positions {
        for k in 0..3 {
            p[k] += amplitude[k] * rng.uniform(-1.0, 1.0);
        }
    }
}

/// A unit sphere obtained by subdividing an icosahedron `level` times and projecting to the
/// sphere. `level = 0` is the icosahedron; there are `10 · 4^level + 2` vertices and
/// `20 · 4^level` faces. Closed, genus 0.
///
/// # Errors
/// `level > 9` (more than 2.6 million vertices) is rejected.
pub fn icosphere(level: u32) -> Result<GeneratedMesh, GeneratorError> {
    if level > 9 {
        return Err(bad(format!("icosphere level must be at most 9, got {level}")));
    }
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut positions: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    for p in &mut positions {
        normalize(p);
    }
    let mut faces: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..level {
        let mut midpoints: HashMap<(u32, u32), u32> = HashMap::with_capacity(faces.len() * 3 / 2);
        let mut next = Vec::with_capacity(faces.len() * 4);
        for &[a, b, c] in &faces {
            let mut mid = |i: u32, j: u32| -> u32 {
                let key = if i < j { (i, j) } else { (j, i) };
                *midpoints.entry(key).or_insert_with(|| {
                    let (p, q) = (positions[i as usize], positions[j as usize]);
                    let mut m = [p[0] + q[0], p[1] + q[1], p[2] + q[2]];
                    normalize(&mut m);
                    positions.push(m);
                    (positions.len() - 1) as u32
                })
            };
            let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
            next.extend_from_slice(&[[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        faces = next;
    }
    Ok(GeneratedMesh { positions, faces })
}

fn normalize(p: &mut [f64; 3]) {
    let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    for x in p {
        *x /= len;
    }
}

/// A torus around the `z` axis with `n_major` segments around the axis and `n_minor` segments
/// around the tube. `n_major · n_minor` vertices, twice as many faces. Closed, genus 1.
///
/// # Errors
/// Each segment count must be at least 3 and `0 < minor_radius < major_radius`.
pub fn torus(
    n_major: usize,
    n_minor: usize,
    major_radius: f64,
    minor_radius: f64,
) -> Result<GeneratedMesh, GeneratorError> {
    if n_major < 3 || n_minor < 3 {
        return Err(bad(format!(
            "torus needs at least 3 segments in each direction, got ({n_major}, {n_minor})"
        )));
    }
    check_count(n_major.checked_mul(n_minor))?;
    check_positive("minor_radius", minor_radius)?;
    check_positive("major_radius", major_radius)?;
    if minor_radius >= major_radius {
        return Err(bad(format!(
            "torus needs minor_radius < major_radius, got {minor_radius} >= {major_radius}"
        )));
    }
    let mut positions = Vec::with_capacity(n_major * n_minor);
    for i in 0..n_major {
        let theta = 2.0 * PI * i as f64 / n_major as f64;
        for j in 0..n_minor {
            let phi = 2.0 * PI * j as f64 / n_minor as f64;
            let ring = major_radius + minor_radius * phi.cos();
            positions.push([ring * theta.cos(), ring * theta.sin(), minor_radius * phi.sin()]);
        }
    }
    let faces = quad_grid_faces(n_major, n_minor, true, true);
    Ok(GeneratedMesh { positions, faces })
}

/// A flat rectangular grid in the plane `z = 0` with `nx × ny` square cells, each split into
/// two triangles. The longer side has length 1. `(nx + 1)(ny + 1)` vertices, vertex `(i, j)` has
/// index `i · (ny + 1) + j` and position `(i h, j h, 0)` with `h = 1 / max(nx, ny)`.
/// One boundary loop, genus 0.
///
/// # Errors
/// `nx` and `ny` must be at least 1.
pub fn grid_disk(nx: usize, ny: usize) -> Result<GeneratedMesh, GeneratorError> {
    bumpy_grid_disk(nx, ny, 0.0)
}

/// Like [`grid_disk`], with height `z = amplitude · sin(3π x) · cos(2π y)` so that the surface
/// has regions of positive and negative curvature.
///
/// # Errors
/// `nx` and `ny` must be at least 1 and `amplitude` finite.
pub fn bumpy_grid_disk(nx: usize, ny: usize, amplitude: f64) -> Result<GeneratedMesh, GeneratorError> {
    if nx < 1 || ny < 1 {
        return Err(bad(format!("grid needs at least 1 cell in each direction, got ({nx}, {ny})")));
    }
    if !amplitude.is_finite() {
        return Err(bad(format!("amplitude must be finite, got {amplitude}")));
    }
    check_count(nx.checked_add(1).and_then(|a| ny.checked_add(1).and_then(|b| a.checked_mul(b))))?;
    let h = 1.0 / nx.max(ny) as f64;
    let mut positions = Vec::with_capacity((nx + 1) * (ny + 1));
    for i in 0..=nx {
        for j in 0..=ny {
            let (x, y) = (i as f64 * h, j as f64 * h);
            let z = if amplitude == 0.0 { 0.0 } else { amplitude * (3.0 * PI * x).sin() * (2.0 * PI * y).cos() };
            positions.push([x, y, z]);
        }
    }
    let faces = quad_grid_faces(nx + 1, ny + 1, false, false);
    Ok(GeneratedMesh { positions, faces })
}

/// A flat annulus in the plane `z = 0` between radii `r_inner` and `r_outer`, with `n_around`
/// segments around and `n_radial` segments across. `n_around · (n_radial + 1)` vertices.
/// Two boundary loops, genus 0, Euler characteristic 0.
///
/// # Errors
/// `n_around >= 3`, `n_radial >= 1` and `0 < r_inner < r_outer` are required.
pub fn annulus(
    n_around: usize,
    n_radial: usize,
    r_inner: f64,
    r_outer: f64,
) -> Result<GeneratedMesh, GeneratorError> {
    check_ring_counts("annulus", n_around, n_radial)?;
    check_positive("r_inner", r_inner)?;
    check_positive("r_outer", r_outer)?;
    if r_inner >= r_outer {
        return Err(bad(format!("annulus needs r_inner < r_outer, got {r_inner} >= {r_outer}")));
    }
    let mut positions = Vec::with_capacity(n_around * (n_radial + 1));
    for i in 0..n_around {
        let theta = 2.0 * PI * i as f64 / n_around as f64;
        for j in 0..=n_radial {
            // Radius decreases with j so that the faces are counter-clockwise seen from +z.
            let r = r_outer + (r_inner - r_outer) * j as f64 / n_radial as f64;
            positions.push([r * theta.cos(), r * theta.sin(), 0.0]);
        }
    }
    let faces = quad_grid_faces(n_around, n_radial + 1, true, false);
    Ok(GeneratedMesh { positions, faces })
}

/// An open cylinder around the `z` axis (no caps) with `n_around` segments around and
/// `n_height` segments along the axis. Same connectivity as [`annulus`]: two boundary loops,
/// genus 0. Faces are counter-clockwise seen from outside.
///
/// # Errors
/// `n_around >= 3`, `n_height >= 1`, and positive `radius` and `height` are required.
pub fn open_cylinder(
    n_around: usize,
    n_height: usize,
    radius: f64,
    height: f64,
) -> Result<GeneratedMesh, GeneratorError> {
    check_ring_counts("cylinder", n_around, n_height)?;
    check_positive("radius", radius)?;
    check_positive("height", height)?;
    let mut positions = Vec::with_capacity(n_around * (n_height + 1));
    for i in 0..n_around {
        let theta = 2.0 * PI * i as f64 / n_around as f64;
        for j in 0..=n_height {
            positions.push([radius * theta.cos(), radius * theta.sin(), height * j as f64 / n_height as f64]);
        }
    }
    let faces = quad_grid_faces(n_around, n_height + 1, true, false);
    Ok(GeneratedMesh { positions, faces })
}

fn check_ring_counts(what: &str, n_around: usize, n_across: usize) -> Result<(), GeneratorError> {
    if n_around < 3 || n_across < 1 {
        return Err(bad(format!(
            "{what} needs at least 3 segments around and 1 across, got ({n_around}, {n_across})"
        )));
    }
    check_count(n_across.checked_add(1).and_then(|a| a.checked_mul(n_around)))
}

fn check_count(n_vertices: Option<usize>) -> Result<(), GeneratorError> {
    match n_vertices {
        Some(n) if n <= (1 << 28) => Ok(()),
        _ => Err(bad("requested mesh is too large (more than 2^28 vertices)".to_string())),
    }
}

/// Triangulates an `ni × nj` lattice of vertices indexed `i · nj + j`. Each cell
/// `(i, j), (i+1, j), (i+1, j+1), (i, j+1)` becomes two triangles sharing the diagonal from
/// `(i, j)` to `(i+1, j+1)`. `wrap_i` / `wrap_j` close the lattice periodically.
fn quad_grid_faces(ni: usize, nj: usize, wrap_i: bool, wrap_j: bool) -> Vec<[u32; 3]> {
    let cells_i = if wrap_i { ni } else { ni - 1 };
    let cells_j = if wrap_j { nj } else { nj - 1 };
    let index = |i: usize, j: usize| ((i % ni) * nj + (j % nj)) as u32;
    let mut faces = Vec::with_capacity(2 * cells_i * cells_j);
    for i in 0..cells_i {
        for j in 0..cells_j {
            let (a, b, c, d) = (index(i, j), index(i + 1, j), index(i + 1, j + 1), index(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    faces
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Signed volume enclosed by a closed mesh; positive for outward orientation.
    fn signed_volume(m: &GeneratedMesh) -> f64 {
        m.faces
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|v| m.positions[v as usize]);
                (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]))
                    / 6.0
            })
            .sum()
    }

    /// Sum of the z components of the face area vectors; positive when faces are
    /// counter-clockwise seen from +z.
    fn area_z(m: &GeneratedMesh) -> f64 {
        m.faces
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|v| m.positions[v as usize]);
                0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
            })
            .sum()
    }

    #[test]
    fn icosphere_counts_and_orientation() {
        for level in 0..4u32 {
            let m = icosphere(level).unwrap();
            let t = m.topology().unwrap();
            assert_eq!(t.n_vertices(), 10 * 4usize.pow(level) + 2);
            assert_eq!(t.n_faces(), 20 * 4usize.pow(level));
            assert_eq!(t.euler_characteristic(), 2);
            assert_eq!(t.genus(), 0);
            assert!(t.boundary_loops().is_empty());
            let vol = signed_volume(&m);
            assert!(vol > 2.5 && vol < 4.0 * PI / 3.0 + 1e-12, "volume {vol}");
            for p in &m.positions {
                let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                assert!((r - 1.0).abs() < 1e-14);
            }
        }
        assert!(icosphere(10).is_err());
    }

    #[test]
    fn torus_topology() {
        let m = torus(12, 7, 2.0, 0.5).unwrap();
        let t = m.topology().unwrap();
        assert_eq!((t.n_vertices(), t.n_faces()), (84, 168));
        assert_eq!(t.euler_characteristic(), 0);
        assert_eq!(t.genus(), 1);
        assert!(t.boundary_loops().is_empty());
        assert!(signed_volume(&m) > 0.0);
        let e = torus(2, 5, 2.0, 0.5).unwrap_err();
        assert_eq!(e.to_string(), "torus needs at least 3 segments in each direction, got (2, 5)");
        assert!(torus(5, 5, 1.0, 1.0).is_err());
        assert!(torus(5, 5, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn grid_disk_topology() {
        for (nx, ny) in [(1, 1), (4, 3), (2, 9)] {
            let m = grid_disk(nx, ny).unwrap();
            let t = m.topology().unwrap();
            assert_eq!(t.n_vertices(), (nx + 1) * (ny + 1));
            assert_eq!(t.n_faces(), 2 * nx * ny);
            assert_eq!(t.euler_characteristic(), 1);
            assert_eq!(t.genus(), 0);
            assert_eq!(t.boundary_loops().len(), 1);
            assert_eq!(t.boundary_loops()[0].len(), 2 * (nx + ny));
            let h = 1.0 / nx.max(ny) as f64;
            assert!((area_z(&m) - (nx as f64 * h) * (ny as f64 * h)).abs() < 1e-12);
        }
        let b = bumpy_grid_disk(6, 6, 0.2).unwrap();
        assert!(b.positions.iter().any(|p| p[2].abs() > 0.05));
        assert_eq!(b.faces, grid_disk(6, 6).unwrap().faces);
        assert!(grid_disk(0, 3).is_err());
    }

    #[test]
    fn annulus_and_cylinder_topology() {
        let a = annulus(9, 3, 0.5, 2.0).unwrap();
        let t = a.topology().unwrap();
        assert_eq!((t.n_vertices(), t.n_faces()), (36, 54));
        assert_eq!(t.euler_characteristic(), 0);
        assert_eq!(t.genus(), 0);
        assert_eq!(t.boundary_loops().len(), 2);
        assert!(area_z(&a) > 0.0);

        let c = open_cylinder(8, 1, 1.0, 2.0).unwrap();
        let t = c.topology().unwrap();
        assert_eq!(t.boundary_loops().len(), 2);
        assert_eq!(t.euler_characteristic(), 0);
        assert_eq!(t.n_boundary_edges(), 16);
        assert!(annulus(2, 1, 0.5, 1.0).is_err());
        assert!(annulus(5, 1, 1.0, 1.0).is_err());
        assert!(open_cylinder(5, 0, 1.0, 1.0).is_err());
    }

    #[test]
    fn perturb_is_deterministic_and_bounded() {
        let base = grid_disk(3, 3).unwrap();
        let mut a = base.positions.clone();
        let mut b = base.positions.clone();
        perturb(&mut a, [0.01, 0.02, 0.0], 7);
        perturb(&mut b, [0.01, 0.02, 0.0], 7);
        assert_eq!(a, b);
        let mut c = base.positions.clone();
        perturb(&mut c, [0.01, 0.02, 0.0], 8);
        assert_ne!(a, c);
        for (p, q) in a.iter().zip(&base.positions) {
            assert!((p[0] - q[0]).abs() <= 0.01 && (p[1] - q[1]).abs() <= 0.02 && p[2] == q[2]);
        }
        // Reference values pin the generator across platforms and releases.
        let mut rng = SplitMix64::new(0);
        assert_eq!(rng.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(rng.next_u64(), 0x6E78_9E6A_A1B9_65F4);
    }
}
