//! Primitives: box, capsule, sphere, torus, punctured torus (Fig. 20), a slab with spheres
//! (Fig. 25), a triangle soup with holes (Fig. 23) and the coarse "H" (Fig. 22).

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};

use super::{Mesh, Rng, bad};
use crate::error::RsError;
use crate::geometry::v3::{norm, sub};

fn from_core(m: ddg_core::generators::GeneratedMesh) -> Mesh {
    Mesh {
        positions: m.positions,
        faces: m.faces,
    }
}

/// The sphere of the given radius around `center`, an icosahedron subdivided `level` times
/// (`20·4^level` faces).
///
/// # Errors
/// `level > 9`, non-positive radius.
pub fn sphere(level: u32, radius: f64, center: [f64; 3]) -> Result<Mesh, RsError> {
    if !(radius.is_finite() && radius > 0.0) {
        return Err(bad("radius", radius, "positive"));
    }
    let m = ddg_core::generators::icosphere(level)
        .map_err(|e| bad("level", level, leak(e.to_string())))?;
    from_core(m).transformed(radius, center).validated()
}

fn leak(s: String) -> &'static str {
    Box::leak(s.into_boxed_str())
}

/// The torus around the `z` axis (`2·n_major·n_minor` faces), genus 1.
///
/// # Errors
/// Segment counts below 3, radii not `0 < minor < major`.
pub fn torus(n_major: usize, n_minor: usize, major: f64, minor: f64) -> Result<Mesh, RsError> {
    let m = ddg_core::generators::torus(n_major, n_minor, major, minor)
        .map_err(|e| bad("torus", format!("{n_major} {n_minor} {major} {minor}"), leak(e.to_string())))?;
    from_core(m).validated()
}

/// The closed box `[-half, half]` with `n` segments per edge (`12 n²` faces), genus 0.
///
/// # Errors
/// `n = 0` or non-positive half extents.
pub fn box_mesh(half: [f64; 3], n: usize) -> Result<Mesh, RsError> {
    if n == 0 || half.iter().any(|&h| !(h.is_finite() && h > 0.0)) {
        return Err(bad("box", format!("{half:?} n={n}"), "n >= 1 and positive half extents"));
    }
    let mut ids: HashMap<[usize; 3], u32> = HashMap::new();
    let mut positions = Vec::new();
    let mut faces = Vec::new();
    let mut id = |c: [usize; 3]| -> u32 {
        *ids.entry(c).or_insert_with(|| {
            positions.push([0, 1, 2].map(|a| -half[a] + 2.0 * half[a] * c[a] as f64 / n as f64));
            (positions.len() - 1) as u32
        })
    };
    for d in 0..3 {
        let (u, v) = ((d + 1) % 3, (d + 2) % 3);
        for s in [0usize, n] {
            for i in 0..n {
                for j in 0..n {
                    let mut c = [[0usize; 3]; 4];
                    for (q, (di, dj)) in [(0, 0), (1, 0), (1, 1), (0, 1)].into_iter().enumerate() {
                        c[q][d] = s;
                        c[q][u] = i + di;
                        c[q][v] = j + dj;
                    }
                    let q = c.map(&mut id);
                    // (u, v, d) is a right-handed frame: u × v = +d
                    if s == n {
                        faces.push([q[0], q[1], q[2]]);
                        faces.push([q[0], q[2], q[3]]);
                    } else {
                        faces.push([q[0], q[2], q[1]]);
                        faces.push([q[0], q[3], q[2]]);
                    }
                }
            }
        }
    }
    Mesh { positions, faces }.validated()
}

/// A capsule along the `z` axis: a cylinder of the given length between two hemispheres.
///
/// # Errors
/// Non-positive sizes, `n_around < 3`, `n_cap < 2`, `n_mid < 1`.
pub fn capsule(
    radius: f64,
    length: f64,
    n_around: usize,
    n_cap: usize,
    n_mid: usize,
) -> Result<Mesh, RsError> {
    if !(radius > 0.0 && length > 0.0 && n_around >= 3 && n_cap >= 2 && n_mid >= 1) {
        return Err(bad("capsule", format!("{radius} {length} {n_around} {n_cap} {n_mid}"), "positive sizes, n_around >= 3, n_cap >= 2, n_mid >= 1"));
    }
    let mut rings: Vec<(f64, f64)> = Vec::new();
    let half = 0.5 * length;
    for m in 1..=n_cap {
        let phi = -FRAC_PI_2 + FRAC_PI_2 * m as f64 / n_cap as f64;
        rings.push((-half + radius * phi.sin(), radius * phi.cos()));
    }
    for q in 1..n_mid {
        rings.push((-half + length * q as f64 / n_mid as f64, radius));
    }
    for m in 0..n_cap {
        let phi = FRAC_PI_2 * m as f64 / n_cap as f64;
        rings.push((half + radius * phi.sin(), radius * phi.cos()));
    }
    let mut positions = vec![[0.0, 0.0, -half - radius]];
    for &(z, rho) in &rings {
        for a in 0..n_around {
            let t = 2.0 * PI * a as f64 / n_around as f64;
            positions.push([rho * t.cos(), rho * t.sin(), z]);
        }
    }
    positions.push([0.0, 0.0, half + radius]);
    let north = (positions.len() - 1) as u32;
    let ring = |r: usize, a: usize| (1 + r * n_around + a % n_around) as u32;
    let mut faces = Vec::new();
    for a in 0..n_around {
        faces.push([0, ring(0, a + 1), ring(0, a)]);
    }
    for r in 0..rings.len() - 1 {
        for a in 0..n_around {
            let (p, q, s, t) = (ring(r, a), ring(r, a + 1), ring(r + 1, a + 1), ring(r + 1, a));
            faces.push([p, q, s]);
            faces.push([p, s, t]);
        }
    }
    let last = rings.len() - 1;
    for a in 0..n_around {
        faces.push([north, ring(last, a), ring(last, a + 1)]);
    }
    Mesh { positions, faces }.validated()
}

/// A torus with a round hole: the faces whose centroid is within `hole_radius` of the outermost
/// point of the equator are removed (Fig. 20). Genus 1, one boundary loop, `χ = −1`.
///
/// # Errors
/// As [`torus`]; a hole that does not leave exactly one boundary loop.
pub fn punctured_torus(
    n_major: usize,
    n_minor: usize,
    major: f64,
    minor: f64,
    hole_radius: f64,
) -> Result<Mesh, RsError> {
    let t = torus(n_major, n_minor, major, minor)?;
    let p0 = [major + minor, 0.0, 0.0];
    let keep: Vec<[u32; 3]> = t
        .faces
        .iter()
        .copied()
        .filter(|f| {
            let c = f.map(|i| t.positions[i as usize]);
            let m = [0, 1, 2].map(|a| (c[0][a] + c[1][a] + c[2][a]) / 3.0);
            norm(sub(m, p0)) > hole_radius
        })
        .collect();
    let m = compact(&t.positions, &keep).validated()?;
    let s = m.surface()?;
    if s.boundary_loops().len() != 1 {
        return Err(bad("hole_radius", hole_radius, "a hole with exactly one boundary loop"));
    }
    Ok(m)
}

/// Keeps only the referenced vertices.
pub fn compact(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> Mesh {
    let mut map = vec![u32::MAX; positions.len()];
    let mut out = Vec::new();
    let mut f2 = Vec::with_capacity(faces.len());
    for f in faces {
        f2.push(f.map(|i| {
            if map[i as usize] == u32::MAX {
                map[i as usize] = out.len() as u32;
                out.push(positions[i as usize]);
            }
            map[i as usize]
        }));
    }
    Mesh {
        positions: out,
        faces: f2,
    }
}

/// Spheres inside a slab: the spheres (separate components, `20·4^level` faces each, radius
/// `radius`) at seeded random, non-overlapping positions at least `gap` from each other and from
/// the walls, and the slab `[-half, half]` as an obstacle (Fig. 25).
#[derive(Debug, Clone)]
pub struct Cobblestones {
    /// The spheres (all components).
    pub spheres: Mesh,
    /// The slab (outward oriented box, obstacle).
    pub slab: Mesh,
    /// The sphere centres.
    pub centers: Vec<[f64; 3]>,
}

/// See [`Cobblestones`].
///
/// # Errors
/// The spheres do not fit (rejection sampling failed) or invalid sizes.
pub fn cobblestones(
    n_spheres: usize,
    radius: f64,
    gap: f64,
    half: [f64; 3],
    level: u32,
    seed: u64,
) -> Result<Cobblestones, RsError> {
    let mut rng = Rng::new(seed);
    let mut centers: Vec<[f64; 3]> = Vec::new();
    let margin = radius + gap;
    for _ in 0..n_spheres {
        let mut placed = false;
        for _ in 0..20000 {
            let c = [0, 1, 2].map(|a| rng.uniform(-half[a] + margin, half[a] - margin));
            if half.iter().any(|&h| h <= margin) {
                break;
            }
            if centers.iter().all(|o| norm(sub(*o, c)) >= 2.0 * radius + gap) {
                centers.push(c);
                placed = true;
                break;
            }
        }
        if !placed {
            return Err(bad("n_spheres", n_spheres, "fewer spheres (they do not fit in the slab)"));
        }
    }
    let mut spheres: Option<Mesh> = None;
    for &c in &centers {
        let s = sphere(level, radius, c)?;
        spheres = Some(match spheres {
            None => s,
            Some(m) => m.merged(&s),
        });
    }
    let spheres = spheres.ok_or_else(|| bad("n_spheres", 0, "at least 1"))?;
    let slab = box_mesh(half, 4)?;
    Ok(Cobblestones {
        spheres,
        slab,
        centers,
    })
}

/// A triangle soup (separate corners, no connectivity) from `mesh` with a fraction of the faces
/// deleted at random (seeded): the polygon soup with holes of Fig. 23.
///
/// # Errors
/// `fraction` outside `[0, 1)`.
pub fn soup_with_holes(
    mesh: &Mesh,
    fraction: f64,
    seed: u64,
) -> Result<Vec<[[f64; 3]; 3]>, RsError> {
    if !(0.0..1.0).contains(&fraction) {
        return Err(bad("fraction", fraction, "in [0, 1)"));
    }
    let mut rng = Rng::new(seed);
    Ok(mesh
        .faces
        .iter()
        .filter(|_| rng.unit() >= fraction)
        .map(|f| f.map(|i| mesh.positions[i as usize]))
        .collect())
}

/// The coarse "H" of Fig. 22: the prism over the 12-gon H (width 3, height 4, depth 1) has 24
/// corner vertices (indices 0..11 on `z = 0`, 12..23 on `z = 1`, same outline order) plus 8 free
/// vertices 24..31 at the centroids of four triangles on each cap: 32 vertices, 60 faces, genus
/// 0. Pinning the 24 corners and penalizing the area deviation (the scenes of
/// `CoarseInterpolation`) lets the 8 free vertices find the minimal-surface-like caps.
pub fn h_shape_coarse() -> Mesh {
    let outline: [[f64; 2]; 12] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.5],
        [2.0, 1.5],
        [2.0, 0.0],
        [3.0, 0.0],
        [3.0, 4.0],
        [2.0, 4.0],
        [2.0, 2.5],
        [1.0, 2.5],
        [1.0, 4.0],
        [0.0, 4.0],
    ];
    let cap = ear_clip(&outline);
    let mut positions: Vec<[f64; 3]> = outline.iter().map(|p| [p[0], p[1], 0.0]).collect();
    positions.extend(outline.iter().map(|p| [p[0], p[1], 1.0]));
    let mut faces: Vec<[u32; 3]> = Vec::new();
    for t in &cap {
        faces.push([t[0], t[2], t[1]]); // bottom, normal -z
        faces.push([t[0] + 12, t[1] + 12, t[2] + 12]); // top, normal +z
    }
    for i in 0..12u32 {
        let j = (i + 1) % 12;
        faces.push([i, j, j + 12]);
        faces.push([i, j + 12, i + 12]);
    }
    // split four triangles of each cap at their centroids
    let mut split = 0;
    let mut out: Vec<[u32; 3]> = Vec::new();
    for f in faces {
        let cap_face = f.iter().all(|&i| i < 12) || f.iter().all(|&i| (12..24).contains(&i));
        if cap_face && split < 8 {
            let c = [0, 1, 2].map(|a| {
                (positions[f[0] as usize][a] + positions[f[1] as usize][a] + positions[f[2] as usize][a]) / 3.0
            });
            positions.push(c);
            let v = (positions.len() - 1) as u32;
            out.push([f[0], f[1], v]);
            out.push([f[1], f[2], v]);
            out.push([f[2], f[0], v]);
            split += 1;
        } else {
            out.push(f);
        }
    }
    Mesh {
        positions,
        faces: out,
    }
}

/// Ear clipping of a simple counter-clockwise polygon.
fn ear_clip(p: &[[f64; 2]]) -> Vec<[u32; 3]> {
    let mut idx: Vec<usize> = (0..p.len()).collect();
    let cross2 = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let mut out = Vec::new();
    while idx.len() > 3 {
        let n = idx.len();
        let mut clipped = false;
        for i in 0..n {
            let (a, b, c) = (idx[(i + n - 1) % n], idx[i], idx[(i + 1) % n]);
            if cross2(p[a], p[b], p[c]) <= 0.0 {
                continue;
            }
            let inside = idx.iter().any(|&q| {
                q != a && q != b && q != c && cross2(p[a], p[b], p[q]) >= 0.0 && cross2(p[b], p[c], p[q]) >= 0.0 && cross2(p[c], p[a], p[q]) >= 0.0
            });
            if !inside {
                out.push([a as u32, b as u32, c as u32]);
                idx.remove(i);
                clipped = true;
                break;
            }
        }
        if !clipped {
            break;
        }
    }
    out.push([idx[0] as u32, idx[1] as u32, idx[2] as u32]);
    out
}

/// Outward-pointing check helper used by tests: the sum of `n_S` weighted by the position
/// relative to the centroid is positive for outward orientation (the signed volume is).
pub fn is_outward(mesh: &Mesh) -> bool {
    mesh.volume() > 0.0
}
