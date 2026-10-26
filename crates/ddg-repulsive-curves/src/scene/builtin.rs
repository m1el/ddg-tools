//! Procedural versions of the applications of `SPEC.md` §14.4, small enough for tests and
//! quick demonstrations. Each returns a [`Scene`]; `size` scales the number of vertices.
//!
//! | name | application | paper |
//! |---|---|---|
//! | [`packing_domain`] | loops packed into a sphere-shaped domain, length growing | Fig. 1, 23 |
//! | [`packing_surface`] | a loop on a torus, length growing | Fig. 24 |
//! | [`packing_boundary`] | planar loops inside a fixed boundary frame | Fig. 22 |
//! | [`graph_drawing`] | a planar graph, straight edges, junctions pinned | Fig. 26 |
//! | [`graph_isometric`] | K5 jittered out of its pentagram drawing, edge lengths kept | Fig. 27 |
//! | [`vascular_tree`] | a branching tree with pinned leaves and kept edge lengths | Fig. 28 |
//! | [`spline_points`] | a loop through control points, growing length | Fig. 30 |
//! | [`spline_tangents`] | an open curve with fixed end tangents and a control point | Fig. 29 |

use std::f64::consts::{PI, TAU};
use std::sync::Arc;

use super::{LengthMode, ObstacleTerm, Scene, SceneSpec, SurfaceVertices};
use crate::error::RcError;
use crate::network::CurveNetwork;
use crate::obstacle::Obstacle;
use crate::penalties::LengthDifference;
use crate::surfaces::{Plane, Torus};

type P = [f64; 3];

/// Names of the generators, for the command line.
pub const NAMES: &[&str] = &[
    "packing-domain",
    "packing-surface",
    "packing-boundary",
    "graph-drawing",
    "graph-isometric",
    "vascular-tree",
    "spline-points",
    "spline-tangents",
];

/// The generator called `name` with `size` (a vertex-count scale, 0 for the default).
///
/// # Errors
/// [`RcError::InvalidOption`] for an unknown name.
pub fn by_name(name: &str, size: usize) -> Result<Scene, RcError> {
    let s = |d: usize| if size == 0 { d } else { size };
    match name {
        "packing-domain" => packing_domain(s(24), 2),
        "packing-surface" => packing_surface(s(48)),
        "packing-boundary" => packing_boundary(s(32), 3),
        "graph-drawing" => graph_drawing(s(10)),
        "graph-isometric" => graph_isometric(s(14)),
        "vascular-tree" => vascular_tree(s(8), 3),
        "spline-points" => spline_points(s(96)),
        "spline-tangents" => spline_tangents(s(60)),
        other => Err(RcError::InvalidOption {
            name: "scene",
            value: other.to_string(),
            requirement: "one of the names in builtin::NAMES",
        }),
    }
}

/// A UV sphere as an obstacle (triangles outward).
pub fn sphere_mesh(center: P, radius: f64, nu: usize, nv: usize) -> Result<Obstacle, RcError> {
    let mut v = vec![[center[0], center[1], center[2] + radius]];
    for i in 1..nv {
        let th = PI * i as f64 / nv as f64;
        for j in 0..nu {
            let ph = TAU * j as f64 / nu as f64;
            v.push([
                center[0] + radius * th.sin() * ph.cos(),
                center[1] + radius * th.sin() * ph.sin(),
                center[2] + radius * th.cos(),
            ]);
        }
    }
    v.push([center[0], center[1], center[2] - radius]);
    let south = (v.len() - 1) as u32;
    let ring = |i: usize, j: usize| (1 + (i - 1) * nu + j % nu) as u32;
    let mut t = Vec::new();
    for j in 0..nu {
        t.push([0, ring(1, j), ring(1, j + 1)]);
        t.push([south, ring(nv - 1, j + 1), ring(nv - 1, j)]);
    }
    for i in 1..nv - 1 {
        for j in 0..nu {
            t.push([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]);
            t.push([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]);
        }
    }
    Obstacle::new(v, t)
}

/// The four walls of the square `[-h, h]^2` in the `xz` plane, extruded to `y in [-h, h]`
/// (a fixed boundary curve for planar packing), `k` cells per side.
pub fn frame_mesh(h: f64, k: usize) -> Result<Obstacle, RcError> {
    let mut v: Vec<P> = Vec::new();
    let mut t = Vec::new();
    let mut wall = |a: [f64; 2], b: [f64; 2]| {
        let base = v.len() as u32;
        for i in 0..=k {
            for j in 0..=k {
                let s = i as f64 / k as f64;
                let y = h * (2.0 * j as f64 / k as f64 - 1.0);
                v.push([a[0] + s * (b[0] - a[0]), y, a[1] + s * (b[1] - a[1])]);
            }
        }
        for i in 0..k {
            for j in 0..k {
                let a0 = base + (i * (k + 1) + j) as u32;
                let b0 = a0 + (k + 1) as u32;
                t.push([a0, b0, b0 + 1]);
                t.push([a0, b0 + 1, a0 + 1]);
            }
        }
    };
    wall([-h, -h], [h, -h]);
    wall([h, -h], [h, h]);
    wall([h, h], [-h, h]);
    wall([-h, h], [-h, -h]);
    Obstacle::new(v, t)
}

fn loop_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect()
}

fn lerp(a: P, b: P, s: f64) -> P {
    [
        a[0] + s * (b[0] - a[0]),
        a[1] + s * (b[1] - a[1]),
        a[2] + s * (b[2] - a[2]),
    ]
}

/// Builds a network from graph vertices (the first `nodes.len()` vertices) and edge paths:
/// each path `(a, b, interior)` becomes the polyline `a, interior..., b`.
pub fn graph_network(
    nodes: &[P],
    paths: &[(usize, usize, Vec<P>)],
) -> Result<(CurveNetwork, Vec<P>), RcError> {
    let mut pos = nodes.to_vec();
    let mut edges = Vec::new();
    for (a, b, interior) in paths {
        let mut prev = *a as u32;
        for p in interior {
            pos.push(*p);
            let id = (pos.len() - 1) as u32;
            edges.push([prev, id]);
            prev = id;
        }
        edges.push([prev, *b as u32]);
    }
    Ok((CurveNetwork::new(pos.len(), edges)?, pos))
}

fn straight(a: P, b: P, k: usize) -> Vec<P> {
    (1..k).map(|i| lerp(a, b, i as f64 / k as f64)).collect()
}

/// Fig. 1, 23: `loops` small circles inside a sphere of radius 1 (a stand-in for the bunny);
/// the length grows by a factor 4 (`rate` 0.04 per step), `n` vertices per loop.
///
/// # Errors
/// Errors of the constructors (none for sensible sizes).
pub fn packing_domain(n: usize, loops: usize) -> Result<Scene, RcError> {
    let mut pos = Vec::new();
    let mut edges = Vec::new();
    for l in 0..loops {
        let ang = TAU * l as f64 / loops as f64;
        let c = [0.35 * ang.cos(), 0.35 * ang.sin(), 0.0];
        let r = 0.18;
        let base = pos.len() as u32;
        for i in 0..n {
            let t = TAU * i as f64 / n as f64;
            // alternate tilts so that the loops are not parallel
            let (ct, st) = (t.cos(), t.sin());
            let q = if l % 2 == 0 {
                [c[0] + r * ct, c[1] + r * st, c[2]]
            } else {
                [c[0] + r * ct, c[1], c[2] + r * st]
            };
            pos.push(q);
        }
        for i in 0..n as u32 {
            edges.push([base + i, base + (i + 1) % n as u32]);
        }
    }
    let net = CurveNetwork::new(pos.len(), edges)?;
    let spec = SceneSpec {
        obstacles: vec![ObstacleTerm {
            mesh: Arc::new(sphere_mesh([0.0; 3], 1.0, 24, 12)?),
            weight: 1.0,
            source: "sphere".into(),
        }],
        barycenter: true,
        length: LengthMode::EdgeLengths { scale: 4.0 },
        schedule_rate: 0.04,
        subdivide_limit: 1,
        iteration_limit: Some(60),
        ..SceneSpec::default()
    };
    Scene::new("packing-domain", net, pos, spec)
}

/// Fig. 24: a loop on the torus `R = 1, r = 0.25` (a (1, 3) curve), edge lengths growing.
///
/// # Errors
/// Errors of the constructors.
pub fn packing_surface(n: usize) -> Result<Scene, RcError> {
    let (rm, rn) = (1.0, 0.25);
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let t = TAU * i as f64 / n as f64;
            let (u, v) = (t, 3.0 * t);
            // the reference's torus has its axis along y
            [
                (rm + rn * v.cos()) * u.cos(),
                rn * v.sin(),
                (rm + rn * v.cos()) * u.sin(),
            ]
        })
        .collect();
    // a (1,3) curve winds three times around the tube: it is long, so start with a small
    // meridian-like loop instead when n is small
    let net = CurveNetwork::new(n, loop_edges(n))?;
    let torus = Torus::new([0.0; 3], rm, rn)?;
    let spec = SceneSpec {
        surface: Some((Arc::new(torus), "torus".into())),
        surface_vertices: SurfaceVertices::All,
        length: LengthMode::EdgeLengths { scale: 1.6 },
        schedule_rate: 0.03,
        iteration_limit: Some(30),
        ..SceneSpec::default()
    };
    Scene::new("packing-surface", net, pos, spec)
}

/// Fig. 22: `k` small circles (`n` vertices each) in the plane `y = 0` inside a square frame
/// of half-width 1; edge lengths grow.
///
/// # Errors
/// Errors of the constructors.
pub fn packing_boundary(n: usize, k: usize) -> Result<Scene, RcError> {
    let mut pos = Vec::new();
    let mut edges = Vec::new();
    for l in 0..k {
        let ang = TAU * l as f64 / k as f64 + 0.3;
        let (cx, cz) = (0.45 * ang.cos(), 0.45 * ang.sin());
        let r = 0.12;
        let base = pos.len() as u32;
        for i in 0..n {
            let t = TAU * i as f64 / n as f64;
            pos.push([cx + r * t.cos(), 0.0, cz + r * t.sin()]);
        }
        for i in 0..n as u32 {
            edges.push([base + i, base + (i + 1) % n as u32]);
        }
    }
    let net = CurveNetwork::new(pos.len(), edges)?;
    let spec = SceneSpec {
        obstacles: vec![ObstacleTerm {
            mesh: Arc::new(frame_mesh(1.0, 6)?),
            weight: 1.0,
            source: "frame".into(),
        }],
        surface: Some((Arc::new(Plane::y_zero()), "yplane".into())),
        surface_vertices: SurfaceVertices::All,
        barycenter: true,
        length: LengthMode::EdgeLengths { scale: 1.5 },
        schedule_rate: 0.02,
        iteration_limit: Some(40),
        ..SceneSpec::default()
    };
    Scene::new("packing-boundary", net, pos, spec)
}

/// Fig. 26: a planar graph in a straight-line layout in the plane `y = 0`: a hexagon with a
/// centre vertex and spokes, every edge bowed so that neighbouring edges crowd each other.
/// Junctions (degree 3 or more) are pinned; the total length is kept; `seg` samples per edge.
///
/// # Errors
/// Errors of the constructors.
pub fn graph_drawing(seg: usize) -> Result<Scene, RcError> {
    let mut nodes: Vec<P> = vec![[0.0; 3]];
    for i in 0..6 {
        let t = TAU * i as f64 / 6.0;
        nodes.push([t.cos(), 0.0, t.sin()]);
    }
    // The junctions sit on a straight-line layout; the edges between them are bowed a little
    // in the plane (taut chains between pinned ends cannot move, and a total-length
    // constraint on straight chains is degenerate).
    let bowed = |a: usize, b: usize, sign: f64| {
        let (p, q) = (nodes[a], nodes[b]);
        let d = [q[0] - p[0], q[2] - p[2]];
        let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
        let nrm = [-d[1] / l, d[0] / l];
        (1..seg)
            .map(|i| {
                let s = i as f64 / seg as f64;
                let w = 0.16 * l * sign * (PI * s).sin();
                [
                    p[0] + s * (q[0] - p[0]) + w * nrm[0],
                    0.0,
                    p[2] + s * (q[2] - p[2]) + w * nrm[1],
                ]
            })
            .collect::<Vec<P>>()
    };
    let mut paths = Vec::new();
    for i in 1..=6usize {
        paths.push((0, i, bowed(0, i, if i % 2 == 0 { 1.0 } else { -1.0 })));
    }
    for i in 1..=6usize {
        let j = i % 6 + 1;
        // the rim bows inwards, towards the spokes
        paths.push((i, j, bowed(i, j, 1.0)));
    }
    let (net, pos) = graph_network(&nodes, &paths)?;
    let spec = SceneSpec {
        // the junctions are pinned
        pins: (0..7).collect(),
        length: LengthMode::TotalLength { scale: 1.0 },
        surface: Some((Arc::new(Plane::y_zero()), "yplane".into())),
        surface_vertices: SurfaceVertices::All,
        iteration_limit: Some(40),
        ..SceneSpec::default()
    };
    Scene::new("graph-drawing", net, pos, spec)
}

/// Fig. 27: K5 on a pentagon, drawn with straight edges (the five diagonals cross in the plane)
/// and lifted out of the plane by a smooth per-edge jitter so that no two edges meet. The
/// junctions are free and every edge keeps its length (an isometric embedding).
///
/// # Errors
/// Errors of the constructors.
pub fn graph_isometric(seg: usize) -> Result<Scene, RcError> {
    let nodes: Vec<P> = (0..5)
        .map(|i| {
            let t = TAU * i as f64 / 5.0 + PI / 2.0;
            [t.cos(), 0.0, t.sin()]
        })
        .collect();
    for seed in 0..64u32 {
        let mut paths = Vec::new();
        let mut ne = 0u32;
        for a in 0..5usize {
            for b in a + 1..5 {
                let amp = 0.08 + 0.05 * (((seed + 3 * ne) % 5) as f64) / 4.0;
                let ph = 0.7 * ((seed * 7 + ne * 13) % 11) as f64;
                let fr = 1.0 + ((seed + ne) % 3) as f64;
                let inter: Vec<P> = (1..seg)
                    .map(|i| {
                        let s = i as f64 / seg as f64;
                        let mut q = lerp(nodes[a], nodes[b], s);
                        q[1] = amp * (PI * s).sin() * (TAU * fr * s + ph).sin();
                        q
                    })
                    .collect();
                paths.push((a, b, inter));
                ne += 1;
            }
        }
        let (net, pos) = graph_network(&nodes, &paths)?;
        let (d, _) = super::min_edge_distance(&net, &pos);
        if d > 0.04 {
            let spec = SceneSpec {
                barycenter: true,
                length: LengthMode::EdgeLengths { scale: 1.0 },
                iteration_limit: Some(40),
                ..SceneSpec::default()
            };
            return Scene::new("graph-isometric", net, pos, spec);
        }
    }
    Err(RcError::InvalidOption {
        name: "graph_isometric",
        value: seg.to_string(),
        requirement: "a jitter with clearance was not found",
    })
}

/// Fig. 28: a binary tree of the given depth in the plane `y = 0` (roots at the bottom),
/// branches bent a little, leaves pinned, junctions free, edge lengths kept, a repulsive
/// obstacle ball in the middle of the canopy (the "organ").
///
/// # Errors
/// Errors of the constructors.
pub fn vascular_tree(seg: usize, depth: usize) -> Result<Scene, RcError> {
    let mut nodes: Vec<P> = vec![[0.0, 0.0, -1.0]];
    let mut paths = Vec::new();
    // (node index, direction angle from +z, length)
    let mut frontier = vec![(0usize, 0.0f64, 0.7f64)];
    for d in 0..depth {
        let mut next = Vec::new();
        for (n, ang, len) in frontier {
            let spread = 0.7 / (1.0 + 0.4 * d as f64);
            if d == 0 {
                // the trunk: one child
                let p = nodes[n];
                let q = [p[0] + len * ang.sin(), 0.0, p[2] + len * ang.cos()];
                nodes.push(q);
                let id = nodes.len() - 1;
                paths.push((n, id, straight(p, q, seg)));
                next.push((id, ang, len * 0.8));
            } else {
                for s in [-1.0, 1.0] {
                    let a = ang + s * spread;
                    let p = nodes[n];
                    let q = [p[0] + len * a.sin(), 0.0, p[2] + len * a.cos()];
                    nodes.push(q);
                    let id = nodes.len() - 1;
                    // bend the polyline slightly sideways
                    let mut inter = straight(p, q, seg);
                    for (i, x) in inter.iter_mut().enumerate() {
                        let w = (PI * (i + 1) as f64 / seg as f64).sin();
                        x[1] = 0.02 * w * s;
                    }
                    paths.push((n, id, inter));
                    next.push((id, a, len * 0.8));
                }
            }
        }
        frontier = next;
    }
    let (net, pos) = graph_network(&nodes, &paths)?;
    let spec = SceneSpec {
        pin_endpoints: true,
        length: LengthMode::EdgeLengths { scale: 1.0 },
        iteration_limit: Some(40),
        obstacles: vec![ObstacleTerm {
            mesh: Arc::new(sphere_mesh([0.3, 0.0, 0.9], 0.05, 8, 6)?),
            weight: 0.05,
            source: "ball".into(),
        }],
        ..SceneSpec::default()
    };
    Scene::new("vascular-tree", net, pos, spec)
}

/// Fig. 30: a closed polygon through 8 control points (pinned) in the plane `y = 0`; the total
/// length grows by a factor 1.4 while `LengthDifference` keeps the spacing even.
///
/// # Errors
/// Errors of the constructors.
pub fn spline_points(n: usize) -> Result<Scene, RcError> {
    let n = n.max(16) / 8 * 8;
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let t = TAU * i as f64 / n as f64;
            let r = 1.0 + 0.35 * (3.0 * t).cos();
            [r * t.cos(), 0.0, r * t.sin()]
        })
        .collect();
    let net = CurveNetwork::new(n, loop_edges(n))?;
    let pins: Vec<u32> = (0..8).map(|k| (k * n / 8) as u32).collect();
    let spec = SceneSpec {
        pins,
        length: LengthMode::TotalLength { scale: 1.4 },
        schedule_rate: 0.02,
        extra_penalties: vec![Arc::new(LengthDifference::new(1.0)?)],
        iteration_limit: Some(40),
        ..SceneSpec::default()
    };
    Scene::new("spline-points", net, pos, spec)
}

/// Fig. 29: an open polyline from `(-1, 0, 0)` to `(1, 0, 0)` with a bend, endpoints pinned
/// with tangents fixed (both pointing into the plane's `+z` direction), one interior control
/// point; all edge lengths grow by 30 % (with only a length penalty or a total-length constraint
/// the end edge collapses and the tangent row degenerates, see the WP13 notes).
///
/// # Errors
/// Errors of the constructors.
pub fn spline_tangents(n: usize) -> Result<Scene, RcError> {
    let n = n.max(12);
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let s = i as f64 / (n - 1) as f64;
            // leaves its start vertically (+z), arrives vertically (-z): a hump
            let x = -1.0 + 2.0 * s;
            let z = 1.2 * (PI * s).sin() * (1.0 + 0.2 * (3.0 * PI * s).sin());
            [x, 0.0, z]
        })
        .collect();
    let edges: Vec<[u32; 2]> = (0..n as u32 - 1).map(|i| [i, i + 1]).collect();
    let net = CurveNetwork::new(n, edges)?;
    let spec = SceneSpec {
        pin_endpoints: true,
        special_tangents: true,
        pins: vec![(n / 2) as u32],
        length: LengthMode::EdgeLengths { scale: 1.3 },
        schedule_rate: 0.02,

        iteration_limit: Some(40),
        ..SceneSpec::default()
    };
    Scene::new("spline-tangents", net, pos, spec)
}
