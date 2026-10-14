//! Simplification of the cut that comes back from an intrinsic triangulation.
//!
//! The pieces of the cut inside a triangle of the caller's triangulation meet at bends: where a
//! segment of the cut crosses an edge of the intrinsic triangulation inside the triangle, or where
//! two segments meet at a point of the export that lies inside it. A bend that turns the cut by
//! little is an artefact of the triangulation the flow ran on. [`Session::merge_bends`] removes
//! such bends (option [`IntrinsicOptions::merge_tolerance`](super::IntrinsicOptions)).

use std::collections::HashMap;

use super::session::Session;

/// A point of the seam-cut mesh.
#[derive(Debug, Clone, Copy)]
pub(super) enum Place {
    Edge { edge: u32, t: f64 },
    Face { face: u32, bary: [f64; 3] },
}

/// A straight piece of the cut inside a face of the seam-cut mesh.
#[derive(Debug, Clone, Copy)]
pub(super) struct Piece {
    pub a: u32,
    pub b: u32,
    pub face: u32,
    /// The regions on the left and on the right of `a -> b`.
    pub left: u32,
    pub right: u32,
    /// How far (a length) the removed bends that this piece stands for lie from it.
    pub slack: f64,
}

/// What the simplification of the cut and the rule for its points did
/// ([`Session::describe_cut_with_report`]).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DescribeReport {
    /// Bends removed ([`IntrinsicOptions::merge_tolerance`](super::IntrinsicOptions)).
    pub n_merged_nodes: usize,
    /// The largest distance, a length in the caller's unit, between the cut as the flow found it
    /// and the cut that was described, over the removed bends (a bound, see the option).
    pub max_displacement: f64,
    /// The same as a fraction of the longest edge of the triangle it happened in.
    pub max_displacement_relative: f64,
    /// Points moved by the rule of [`IntrinsicOptions::min_offset`](super::IntrinsicOptions)
    /// (also in `n_clamped_crossings`, with those the flow moved).
    pub n_moved_points: usize,
}

impl Session {
    /// Removes bends inside triangles of the seam-cut mesh. A bend `B` between the pieces `A-B`
    /// and `B-C` of one triangle is removed (the pieces replaced by `A-C`) when
    ///
    /// * the points `A` and `C` do not lie on one side of the triangle (the new piece would
    ///   run along it);
    /// * it is a point inside the triangle with exactly these two pieces (never a junction of
    ///   three or more), and both pieces have the same regions on their sides;
    /// * the distance of `B` from `A-C`, added to the slack the two pieces already carry (the
    ///   bends removed earlier that they stand for), is at most `tolerance` times the longest
    ///   edge of the triangle: the cut is nowhere moved further than that;
    /// * no other point of the triangle lies in the closed triangle `A B C` and no other piece
    ///   crosses `A-C`. Since `A B C` lies inside the triangle and its corners are extreme
    ///   points, no vertex of the mesh is inside it: the cut keeps passing on the same side of
    ///   every vertex.
    pub(super) fn merge_bends(
        &self,
        places: &[Place],
        pieces: &mut Vec<Piece>,
        tolerance: f64,
        report: &mut DescribeReport,
    ) {
        let x = &self.flow_positions;
        let faces = self.flow_mesh.faces();
        // The barycentric coordinates of a place in a face, if it lies on the closed face.
        let bary_in = |place: &Place, face: u32| -> Option<[f64; 3]> {
            let f = faces[face as usize];
            match *place {
                Place::Face { face: g, bary } => (g == face).then_some(bary),
                Place::Edge { edge, t } => {
                    let [a, b] = self.flow_mesh.edge_vertices(edge as usize);
                    let mut out = [0.0; 3];
                    let mut hit = 0;
                    for k in 0..3 {
                        if f[k] as usize == a {
                            out[k] = 1.0 - t;
                            hit += 1;
                        } else if f[k] as usize == b {
                            out[k] = t;
                            hit += 1;
                        }
                    }
                    (hit == 2).then_some(out)
                }
            }
        };
        let point3 = |face: u32, b: [f64; 3]| -> [f64; 3] {
            let f = faces[face as usize];
            [0, 1, 2].map(|d| (0..3).map(|k| b[k] * x[f[k] as usize][d]).sum())
        };
        let cross = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        };
        let xy = |b: [f64; 3]| [b[1], b[2]];
        let mut alive = vec![true; pieces.len()];
        for _pass in 0..200 {
            let mut changed = false;
            // The pieces at every point inside a face.
            let mut at: HashMap<u32, Vec<usize>> = HashMap::new();
            let mut touched = vec![false; places.len()];
            for (i, p) in pieces.iter().enumerate() {
                if !alive[i] {
                    continue;
                }
                for g in [p.a, p.b] {
                    touched[g as usize] = true;
                    if matches!(places[g as usize], Place::Face { .. }) {
                        at.entry(g).or_default().push(i);
                    }
                }
            }
            let mut nodes: Vec<u32> = at.keys().copied().collect();
            nodes.sort_unstable();
            for b_node in nodes {
                let list = &at[&b_node];
                if list.len() != 2 || list.iter().any(|&i| !alive[i]) {
                    continue;
                }
                let (i, j) = (list[0], list[1]);
                let (pi, pj) = (pieces[i], pieces[j]);
                let a_node = if pi.a == b_node { pi.b } else { pi.a };
                let c_node = if pj.a == b_node { pj.b } else { pj.a };
                if a_node == c_node || pi.face != pj.face {
                    continue;
                }
                // Walking A -> B -> C the regions on the left and on the right must agree.
                let lr = |p: &Piece, forward: bool| {
                    if forward {
                        (p.left, p.right)
                    } else {
                        (p.right, p.left)
                    }
                };
                let (li, lj) = (lr(&pi, pi.a == a_node), lr(&pj, pj.a == b_node));
                if li != lj {
                    continue;
                }
                let face = pi.face;
                let (Some(ba), Some(bb), Some(bc)) = (
                    bary_in(&places[a_node as usize], face),
                    bary_in(&places[b_node as usize], face),
                    bary_in(&places[c_node as usize], face),
                ) else {
                    continue;
                };
                // A shortcut between two points of one side would run along that side.
                if (0..3).any(|k| ba[k] < 1e-12 && bc[k] < 1e-12) {
                    continue;
                }
                // Distance of B from the segment A C in space.
                let (pa, pb, pc) = (point3(face, ba), point3(face, bb), point3(face, bc));
                let ac = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
                let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                let len2: f64 = ac.iter().map(|v| v * v).sum();
                if len2 == 0.0 {
                    continue;
                }
                let s =
                    (ab.iter().zip(&ac).map(|(p, q)| p * q).sum::<f64>() / len2).clamp(0.0, 1.0);
                let d = (0..3)
                    .map(|k| (ab[k] - s * ac[k]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                let f = faces[face as usize];
                let longest = (0..3)
                    .map(|k| {
                        let (u, v) = (x[f[k] as usize], x[f[(k + 1) % 3] as usize]);
                        (0..3).map(|c| (u[c] - v[c]).powi(2)).sum::<f64>().sqrt()
                    })
                    .fold(0.0, f64::max);
                let slack = pi.slack.max(pj.slack) + d;
                if slack > tolerance * longest {
                    continue;
                }
                // In the plane of the face: no other point inside the triangle A B C, no piece
                // crossing A C.
                let (qa, qb, qc) = (xy(ba), xy(bb), xy(bc));
                let sign = cross(qa, qb, qc).signum();
                let mut blocked = false;
                for (g, place) in places.iter().enumerate() {
                    let g = g as u32;
                    if !touched[g as usize] || [a_node, b_node, c_node].contains(&g) {
                        continue;
                    }
                    let Some(bg) = bary_in(place, face) else {
                        continue;
                    };
                    let q = xy(bg);
                    let eps = 1e-12;
                    if cross(qa, qb, q) * sign >= -eps
                        && cross(qb, qc, q) * sign >= -eps
                        && cross(qc, qa, q) * sign >= -eps
                    {
                        blocked = true;
                        break;
                    }
                }
                if !blocked {
                    for (k, other) in pieces.iter().enumerate() {
                        if !alive[k] || k == i || k == j || other.face != face {
                            continue;
                        }
                        if [other.a, other.b].contains(&a_node)
                            || [other.a, other.b].contains(&c_node)
                        {
                            continue;
                        }
                        let (Some(u), Some(v)) = (
                            bary_in(&places[other.a as usize], face),
                            bary_in(&places[other.b as usize], face),
                        ) else {
                            continue;
                        };
                        let (u, v) = (xy(u), xy(v));
                        let (d1, d2) = (cross(qa, qc, u), cross(qa, qc, v));
                        let (d3, d4) = (cross(u, v, qa), cross(u, v, qc));
                        if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                            blocked = true;
                            break;
                        }
                    }
                }
                if blocked {
                    continue;
                }
                let (left, right) = li;
                alive[i] = false;
                alive[j] = false;
                pieces.push(Piece {
                    a: a_node,
                    b: c_node,
                    face,
                    left,
                    right,
                    slack,
                });
                alive.push(true);
                report.n_merged_nodes += 1;
                report.max_displacement = report.max_displacement.max(slack);
                report.max_displacement_relative =
                    report.max_displacement_relative.max(slack / longest);
                changed = true;
            }
            if !changed {
                break;
            }
        }
        let kept: Vec<Piece> = pieces
            .iter()
            .zip(&alive)
            .filter(|(_, a)| **a)
            .map(|(p, _)| *p)
            .collect();
        *pieces = kept;
    }
}
