//! Shapes meshed by marching tetrahedra: genus-g surfaces (Fig. 15), tube surfaces around curve
//! networks (handcuffs, Fig. 12 and 16) and knotted tunnels (Fig. 12). The Fig. 12 meshes of the
//! paper are not published; these are analogues with the paper's face counts (`SPEC.md` D31).

use std::f64::consts::TAU;

use super::check::clearance;
use super::field::{self, Field};
use super::tetra::mesh_with_faces;
use super::{Mesh, bad};
use crate::error::RsError;
use crate::geometry::v3::{dot, norm, scale};

type P = [f64; 3];

/// The face counts of the five Fig. 12 meshes: linked, doubly-wound and 6×-wound handcuffs,
/// figure-8 tunnel, trefoil tunnel.
pub const FIG12_FACES: [usize; 5] = [25_098, 26_306, 21_712, 45_198, 23_206];

fn expect_topology(mesh: &Mesh, chi: i64, what: &'static str) -> Result<(), RsError> {
    let s = mesh.surface()?;
    if s.n_components() != 1 || s.euler_characteristic() != chi {
        return Err(bad(
            what,
            format!(
                "{} faces: {} components, chi {} (expected 1, {chi})",
                mesh.faces.len(),
                s.n_components(),
                s.euler_characteristic()
            ),
            "more faces (the tube is not resolved by the grid)",
        ));
    }
    Ok(())
}

/// A surface with the core curves that generated it.
#[derive(Debug, Clone)]
pub struct TubeShape {
    /// The surface.
    pub mesh: Mesh,
    /// The closed core curves (polylines), e.g. the two loops of the handcuffs.
    pub cores: Vec<Vec<P>>,
    /// The tube radius.
    pub radius: f64,
}

/// Genus-`g` surface with the handles in a row (Fig. 15): a rounded slab (half extents
/// `(1.2(g−1) + 1.2, 1.2, 0.5)`) with `g` round holes of radius 0.5 spaced 2.4 apart; `g = 0` is
/// the slab without holes.
/// About `faces` triangles.
///
/// # Errors
/// `faces < 100`, `g > 40`, or meshing errors.
pub fn genus_surface(g: u32, faces: usize) -> Result<Mesh, RsError> {
    if g > 40 {
        return Err(bad("g", g, "at most 40"));
    }
    let (sp, rh, t) = (2.4, 0.5, 0.5);
    let hx = if g == 0 { 1.2 } else { (g as f64 - 1.0) * sp / 2.0 + rh + 0.7 };
    let hy = if g == 0 { 1.0 } else { rh + 0.7 };
    let mut f: Field = field::rounded_box([0.0; 3], [hx, hy, t], 0.45);
    for i in 0..g {
        let x = (i as f64 - (g as f64 - 1.0) / 2.0) * sp;
        f = field::difference(f, field::cylinder([x, 0.0, 0.0], [0.0, 0.0, 1.0], rh), 0.15);
    }
    let m = 0.8;
    mesh_with_faces(&f, [-hx - m, -hy - m, -t - m], [hx + m, hy + m, t + m], faces)
}

/// The two parallel `(1, k)` curves of the torus link `T(2, 2k)` on the torus of radii
/// `major, minor`, `n` points each: their linking number is `k` (the Hopf link for `k = 1`).
pub fn torus_link_cores(k: u32, major: f64, minor: f64, n: usize) -> [Vec<P>; 2] {
    let curve = |phase: f64| -> Vec<P> {
        (0..n)
            .map(|i| {
                let t = TAU * i as f64 / n as f64;
                let th = k as f64 * t + phase;
                let ring = major + minor * th.cos();
                [ring * t.cos(), ring * t.sin(), minor * th.sin()]
            })
            .collect()
    };
    [curve(0.0), curve(std::f64::consts::PI)]
}

/// A tube of radius `radius` around the given curves `(points, closed)` (each blended into the
/// others with width `radius`), meshed with about `faces` triangles. The curves must form a
/// graph whose tube is a manifold surface (the caller makes sure the tube does not touch
/// itself).
///
/// # Errors
/// Meshing errors; `radius ≤ 0`.
pub fn graph_tube(curves: &[(Vec<P>, bool)], radius: f64, faces: usize) -> Result<Mesh, RsError> {
    if !(radius.is_finite() && radius > 0.0) {
        return Err(bad("radius", radius, "positive"));
    }
    let f = field::tube_network(
        curves.iter().map(|(p, c)| (p.clone(), *c, radius)).collect(),
        radius,
    );
    let (lo, hi) = bounds(curves.iter().flat_map(|(p, _)| p.iter().copied()));
    let ext = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max);
    let m = 2.0 * radius + ext / 25.0;
    mesh_with_faces(&f, lo.map(|x| x - m), hi.map(|x| x + m), faces)
}

fn bounds(pts: impl Iterator<Item = P>) -> (P, P) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

/// Handcuffs of winding number `k ≥ 1`: a tube (genus 2) around two unknotted loops, the cores of
/// the torus link `T(2, 2k)` (linking number `k`), joined by a short bridge across the torus
/// tube. `k = 1`: linked handcuffs; `k = 2`: doubly wound; `k = 6`: 6×-wound (the Fig. 12 meshes
/// are not published; these are analogues, D31). About `faces` triangles.
///
/// # Errors
/// `k = 0`, `k > 12`, meshing errors; a mesh whose topology is not the designed one (genus 2,
/// one component) because the grid does not resolve the tube: ask for more faces (about 20 000
/// for `k = 6`).
pub fn handcuffs(k: u32, faces: usize) -> Result<TubeShape, RsError> {
    if k == 0 || k > 12 {
        return Err(bad("k", k, "between 1 and 12"));
    }
    let (major, minor, n) = (3.0, 1.0, 360);
    let [a, b] = torus_link_cores(k, major, minor, n);
    let excl = (n / (4 * k as usize)).max(2);
    let gap = clearance(&a, &b, false, 0)
        .min(clearance(&a, &a, true, excl))
        .min(clearance(&b, &b, true, excl));
    let radius = 0.3 * gap;
    // bridge across the torus tube at t = 0 (angle 0 to angle π)
    let bridge = vec![a[0], b[0]];
    let curves = vec![(a.clone(), true), (b.clone(), true), (bridge, false)];
    let mesh = graph_tube(&curves, radius, faces)?;
    expect_topology(&mesh, -2, "handcuffs")?;
    Ok(TubeShape {
        mesh,
        cores: vec![a, b],
        radius,
    })
}

/// Which knot a tunnel follows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Knot {
    /// The trefoil `(sin t + 2 sin 2t, cos t − 2 cos 2t, −sin 3t)`.
    Trefoil,
    /// The figure-eight knot `((2 + cos 2t) cos 3t, (2 + cos 2t) sin 3t, sin 4t)`.
    FigureEight,
}

/// The knot curve with `n` points.
pub fn knot_curve(knot: Knot, n: usize) -> Vec<P> {
    (0..n)
        .map(|i| {
            let t = TAU * i as f64 / n as f64;
            match knot {
                Knot::Trefoil => [t.sin() + 2.0 * (2.0 * t).sin(), t.cos() - 2.0 * (2.0 * t).cos(), -(3.0 * t).sin()],
                Knot::FigureEight => {
                    let r = 2.0 + (2.0 * t).cos();
                    [r * (3.0 * t).cos(), r * (3.0 * t).sin(), (4.0 * t).sin()]
                }
            }
        })
        .collect()
}

/// A knotted tunnel (Fig. 12): the solid `ball ∩ {x·d ≤ c}` (a flat-topped ball) from which a
/// tube along `knot` is removed. The direction `d` is chosen (deterministically, from 400
/// directions) so that exactly one lobe of the knot sticks out of the flat top: the knot enters
/// and leaves the body through the top, and the surface is a knotted torus (genus 1, one
/// component, `χ = 0`). `cores` holds the knot. About `faces` triangles.
///
/// # Errors
/// No direction separates one lobe by at least 4 tube radii; meshing errors.
pub fn knotted_tunnel(knot: Knot, faces: usize) -> Result<TubeShape, RsError> {
    let n = 600;
    let pts = knot_curve(knot, n);
    let gap = clearance(&pts, &pts, true, n / 16);
    let radius = 0.3 * gap;
    // direction with the largest height gap between the two highest lobes
    let w = n / 20;
    let mut best = (f64::NEG_INFINITY, [0.0, 0.0, 1.0], 0.0);
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    for s in 0..400 {
        let z = 1.0 - 2.0 * (s as f64 + 0.5) / 400.0;
        let r = (1.0 - z * z).sqrt();
        let d = [r * (golden * s as f64).cos(), r * (golden * s as f64).sin(), z];
        let h: Vec<f64> = pts.iter().map(|p| dot(*p, d)).collect();
        let mut maxima: Vec<f64> = (0..n)
            .filter(|&i| (1..=w).all(|o| h[i] >= h[(i + o) % n] && h[i] >= h[(i + n - o) % n]))
            .map(|i| h[i])
            .collect();
        maxima.sort_by(|a, b| b.total_cmp(a));
        // lobes: collapse maxima that belong to the same hump (within a tube radius in height
        // and adjacent) by taking distinct ones separated by > 0.5 radius
        maxima.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
        if maxima.len() >= 2 {
            let gap_h = maxima[0] - maxima[1];
            if gap_h > best.0 {
                best = (gap_h, d, 0.5 * (maxima[0] + maxima[1]));
            }
        }
    }
    let (gap_h, d, plane) = best;
    if gap_h < 4.0 * radius {
        return Err(bad("knot", gap_h / radius, "a lobe separated by at least 4 tube radii"));
    }
    let rs = pts.iter().map(|p| norm(*p)).fold(0.0, f64::max) + 3.0 * radius;
    let ball = field::sphere([0.0; 3], rs);
    let cut = field::half_space(scale(d, plane), d);
    let body = field::intersection(ball, cut, 0.0);
    let tubef = field::tube(pts.clone(), true, radius);
    let f = field::difference(body, tubef, radius * 0.5);
    let m = radius + rs * 0.1;
    let mesh = mesh_with_faces(&f, [-rs - m; 3], [rs + m; 3], faces)?;
    expect_topology(&mesh, 0, "tunnel")?;
    Ok(TubeShape {
        mesh,
        cores: vec![pts],
        radius,
    })
}

/// The five Fig. 12 analogues at the paper's face counts: linked, doubly-wound, 6×-wound
/// handcuffs, figure-8 tunnel, trefoil tunnel.
///
/// # Errors
/// As the constructors.
pub fn fig12(index: usize) -> Result<TubeShape, RsError> {
    let faces = *FIG12_FACES.get(index).ok_or_else(|| bad("index", index, "0..5"))?;
    match index {
        0 => handcuffs(1, faces),
        1 => handcuffs(2, faces),
        2 => handcuffs(6, faces),
        3 => knotted_tunnel(Knot::FigureEight, faces),
        _ => knotted_tunnel(Knot::Trefoil, faces),
    }
}
