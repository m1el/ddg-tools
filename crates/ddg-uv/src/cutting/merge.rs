//! Simplification of the cut that comes back from an intrinsic triangulation.
//!
//! The pieces of the cut inside a triangle of the caller's triangulation meet at bends: where a
//! segment of the cut crosses an edge of the intrinsic triangulation inside the triangle, or where
//! two segments meet at a point of the export that lies inside it. A bend that turns the cut by
//! little is an artefact of the triangulation the flow ran on. [`Session::merge_bends`] removes
//! such bends (option [`IntrinsicOptions::merge_tolerance`](super::IntrinsicOptions)), and the
//! excursions of the cut across a side of a triangle whose points the rule of
//! [`IntrinsicOptions::min_offset`](super::IntrinsicOptions) moved.
//!
//! Every point has two places: where the way back put it (`found`) and where the rule of
//! `min_offset` moved it (`places`, what is described and inserted). The displacement of a merge
//! is measured at the found places, so that it is the displacement of the cut the flow found;
//! the tests that a merge crosses nothing are made at the moved places, which are inserted.

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
    /// and the cut that was described, over the removed bends (a bound, see the option). It is
    /// measured before the rule of `min_offset` moves the points, which moves them by its own
    /// amount on top of this (counted in `n_moved_points`).
    pub max_displacement: f64,
    /// The same as a fraction of the longest edge of the triangle it happened in.
    pub max_displacement_relative: f64,
    /// Points moved by the rule of [`IntrinsicOptions::min_offset`](super::IntrinsicOptions)
    /// (also in `n_clamped_crossings`, with those the flow moved).
    pub n_moved_points: usize,
    /// Of the bends removed, those that were excursions of the cut across a side of a triangle
    /// (the cut crosses the side, bends inside the triangle and crosses the same side back):
    /// they are removed with their two crossings of the side, and the cut stays in the
    /// neighbouring triangle (see [`Session::merge_bends`]).
    pub n_removed_excursions: usize,
}

/// A bend `B` whose two pieces `A-B`, `B-C` end on one side of the triangle `face`.
struct Excursion {
    b: u32,
    a: u32,
    c: u32,
    /// The pieces `A-B` and `B-C`.
    pieces: [usize; 2],
    face: u32,
    /// The regions on the left and on the right walking `A -> B -> C`.
    regions: (u32, u32),
}

/// What [`Session::remove_excursion`] replaces.
struct Shortcut {
    /// The pieces from `X` to `Y`.
    dead: Vec<usize>,
    /// The piece `X-Y`.
    piece: Piece,
    /// The edge length the displacement is measured against.
    longest: f64,
}

impl Session {
    /// Removes bends inside triangles of the seam-cut mesh. A bend `B` between the pieces `A-B`
    /// and `B-C` of one triangle is removed (the pieces replaced by `A-C`) when
    ///
    /// * the points `A` and `C` do not lie on one side of the triangle (the new piece would
    ///   run along it; such a bend is an excursion across that side, removed by
    ///   [`Session::remove_excursion`] if the rule of `min_offset` moved `A`, `B` or `C`);
    /// * it is a point inside the triangle with exactly these two pieces (never a junction of
    ///   three or more), and both pieces have the same regions on their sides;
    /// * the distance of `B` from `A-C`, added to the slack the two pieces already carry (the
    ///   bends removed earlier that they stand for), is at most `tolerance` times the longest
    ///   edge of the triangle: the cut is nowhere moved further than that. The distance is taken
    ///   at the places where the way back put the three points (`found`); at the places the rule
    ///   of `min_offset` moved them to, a bend pushed off a side would count the rule's motion
    ///   as the merge's and stay (with `min_offset` 0.2 and the optimal Delaunay triangulation
    ///   this kept about 6 % of the bends of the spot atlas, and the unwrap inverted parts at
    ///   them);
    /// * no other point of the triangle lies in the closed triangle `A B C` and no other piece
    ///   crosses `A-C` (at the moved places, which are inserted). Since `A B C` lies inside the triangle and its corners are extreme
    ///   points, no vertex of the mesh is inside it: the cut keeps passing on the same side of
    ///   every vertex.
    pub(super) fn merge_bends(
        &self,
        places: &[Place],
        found: &[Place],
        pieces: &mut Vec<Piece>,
        tolerance: f64,
        report: &mut DescribeReport,
    ) {
        let x = &self.flow_positions;
        let faces = self.flow_mesh.faces();
        let bary_in = |place: &Place, face: u32| self.bary_in(place, face);
        let point3 = |face: u32, b: [f64; 3]| self.point3(face, b);
        let cross = cross2;
        let xy = |b: [f64; 3]| [b[1], b[2]];
        let mut alive = vec![true; pieces.len()];
        for _pass in 0..200 {
            let mut changed = false;
            // The pieces at every point; the bends are points inside a face.
            let mut at: HashMap<u32, Vec<usize>> = HashMap::new();
            let mut touched = vec![false; places.len()];
            for (i, p) in pieces.iter().enumerate() {
                if !alive[i] {
                    continue;
                }
                for g in [p.a, p.b] {
                    touched[g as usize] = true;
                    at.entry(g).or_default().push(i);
                }
            }
            let mut nodes: Vec<u32> = at
                .keys()
                .copied()
                .filter(|&g| matches!(places[g as usize], Place::Face { .. }))
                .collect();
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
                // A shortcut between two points of one side would run along that side: the
                // bend is an excursion of the cut across that side, removed (if at all) from the
                // other face.
                // Only an excursion whose bend or crossings the rule of `min_offset` moved is
                // removed: that rule is what turns it into a spike (a bend close to the side) or
                // makes its two crossings one point (both close to a vertex).
                if (0..3).any(|k| ba[k] < 1e-12 && bc[k] < 1e-12) {
                    if ![a_node, b_node, c_node]
                        .iter()
                        .any(|&g| moved_by_the_rule(&places[g as usize], &found[g as usize]))
                    {
                        continue;
                    }
                    let excursion = Excursion {
                        b: b_node,
                        a: a_node,
                        c: c_node,
                        pieces: [i, j],
                        face,
                        regions: li,
                    };
                    if let Some(merged) = self.remove_excursion(
                        &excursion, places, found, pieces, &alive, &at, &touched, tolerance,
                    ) {
                        for k in merged.dead {
                            alive[k] = false;
                        }
                        let longest = merged.longest;
                        let slack = merged.piece.slack;
                        pieces.push(merged.piece);
                        alive.push(true);
                        report.n_merged_nodes += 1;
                        report.n_removed_excursions += 1;
                        report.max_displacement = report.max_displacement.max(slack);
                        report.max_displacement_relative =
                            report.max_displacement_relative.max(slack / longest);
                        changed = true;
                    }
                    continue;
                }
                // Distance of B from the segment A C in space, at the places where the way back
                // put the three points (before the rule of `min_offset` moved them): the
                // displacement is that of the cut the flow found. Measured at the moved places,
                // a bend that the rule pushed away from a side would count the rule's motion as
                // the merge's, and stay.
                let (Some(fa), Some(fb), Some(fc)) = (
                    bary_in(&found[a_node as usize], face),
                    bary_in(&found[b_node as usize], face),
                    bary_in(&found[c_node as usize], face),
                ) else {
                    continue;
                };
                let (pa, pb, pc) = (point3(face, fa), point3(face, fb), point3(face, fc));
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

    /// Removes the excursion `e` of the cut across a side `s` of its triangle `T`: the cut comes
    /// from `X` in the neighbouring triangle `T'`, crosses `s` at `A`, bends at `B` inside `T`
    /// and crosses `s` back at `C` to `Y` in `T'`. If `X` and `Y` lie on one side of `T'` too,
    /// the excursion continues across that side, and the walk goes on outwards in the same way.
    /// The pieces from the outermost `X` to the outermost `Y` become one piece `X-Y` of the
    /// outermost triangle when
    ///
    /// * every pair of crossings (`A` and `C`, then the `X` and `Y` of each step inwards) lies
    ///   on one edge, has exactly two pieces at each point, and the next pieces lie in the other
    ///   face of that edge; the outermost `X ≠ Y` are not on one side of their triangle; all
    ///   pieces have the same regions on their sides; at most 8 steps;
    /// * the largest distance of the removed points from `X-Y` (at the places where the way
    ///   back put them), added to the slack the removed pieces carry, is at most `tolerance`
    ///   times the longest edge of the smallest of the triangles;
    /// * no other point of the cut lies in the closed triangle `A B C` of `T` or in the convex
    ///   hull of the four points of the walk in each further triangle, and no other piece of the
    ///   outermost triangle crosses `X-Y`. These regions lie in their triangles and have no
    ///   vertex of the mesh inside, so the cut keeps passing on the same side of every vertex.
    ///
    /// Without this an excursion stays: the shortcut `A-C` would run along the side. When the
    /// cut only touches the triangle its bend is very close to the side, and the rule of
    /// `min_offset` moves it into the triangle by that fraction of the height: the excursion
    /// becomes a spike whose parts the flattening of the chart inverts.
    #[allow(clippy::too_many_arguments)]
    fn remove_excursion(
        &self,
        e: &Excursion,
        places: &[Place],
        found: &[Place],
        pieces: &[Piece],
        alive: &[bool],
        at: &HashMap<u32, Vec<usize>>,
        touched: &[bool],
        tolerance: f64,
    ) -> Option<Shortcut> {
        const MAX_STEPS: usize = 8;
        // The other piece at a crossing with exactly two.
        let other = |node: u32, mine: usize| -> Option<usize> {
            let list = at.get(&node)?;
            if list.len() != 2 || list.iter().any(|&k| !alive[k]) {
                return None;
            }
            list.iter().copied().find(|&k| k != mine)
        };
        let lr = |p: &Piece, forward: bool| {
            if forward {
                (p.left, p.right)
            } else {
                (p.right, p.left)
            }
        };
        let plane = |b: [f64; 3]| [b[1], b[2]];
        let mut swept: Vec<(u32, Vec<[f64; 2]>)> = vec![(
            e.face,
            vec![
                plane(self.bary_in(&places[e.a as usize], e.face)?),
                plane(self.bary_in(&places[e.b as usize], e.face)?),
                plane(self.bary_in(&places[e.c as usize], e.face)?),
            ],
        )];
        let mut dead = vec![e.pieces[0], e.pieces[1]];
        let mut removed = vec![e.a, e.b, e.c];
        let (mut a, mut c) = (e.a, e.c);
        let (mut piece_a, mut piece_c) = (e.pieces[0], e.pieces[1]);
        let mut inner = e.face;
        let (x, y, far) = loop {
            let (Place::Edge { edge: ea, .. }, Place::Edge { edge: ec, .. }) =
                (places[a as usize], places[c as usize])
            else {
                return None;
            };
            if ea != ec || swept.len() > MAX_STEPS {
                return None;
            }
            let (ia, ic) = (other(a, piece_a)?, other(c, piece_c)?);
            let (xa, cy) = (pieces[ia], pieces[ic]);
            let far = xa.face;
            if cy.face != far || far == inner {
                return None;
            }
            let x = if xa.a == a { xa.b } else { xa.a };
            let y = if cy.a == c { cy.b } else { cy.a };
            if x == y || x == c || y == a {
                return None;
            }
            // Walking X -> A -> ... -> C -> Y the regions on the left and on the right agree.
            if lr(&xa, xa.a == x) != e.regions || lr(&cy, cy.a == c) != e.regions {
                return None;
            }
            dead.extend([ia, ic]);
            let (bx, by) = (
                self.bary_in(&places[x as usize], far)?,
                self.bary_in(&places[y as usize], far)?,
            );
            swept.push((
                far,
                vec![
                    plane(bx),
                    plane(self.bary_in(&places[a as usize], far)?),
                    plane(self.bary_in(&places[c as usize], far)?),
                    plane(by),
                ],
            ));
            if !(0..3).any(|k| bx[k] < 1e-12 && by[k] < 1e-12) {
                break (x, y, far);
            }
            // X and Y on one side of this triangle: the excursion continues across it.
            removed.extend([x, y]);
            (a, c, piece_a, piece_c, inner) = (x, y, ia, ic, far);
        };
        // The displacement, at the places where the way back put the points.
        let (fx, fy) = (
            self.position(&found[x as usize]),
            self.position(&found[y as usize]),
        );
        let d = removed
            .iter()
            .map(|&g| distance_to_segment(self.position(&found[g as usize]), fx, fy))
            .fold(0.0, f64::max);
        let carried = dead.iter().map(|&k| pieces[k].slack).fold(0.0, f64::max);
        let slack = carried + d;
        let longest = swept
            .iter()
            .map(|(face, _)| self.longest_edge(*face))
            .fold(f64::INFINITY, f64::min);
        if slack > tolerance * longest || slack.is_nan() {
            return None;
        }
        // No other point in the swept regions.
        let mut own = removed.clone();
        own.extend([x, y]);
        for (g, place) in places.iter().enumerate() {
            if !touched[g] || own.contains(&(g as u32)) {
                continue;
            }
            for (face, region) in &swept {
                if let Some(b) = self.bary_in(place, *face)
                    && in_hull(region, plane(b))
                {
                    return None;
                }
            }
        }
        // No other piece of the outermost triangle crosses X-Y.
        let (qx, qy) = (
            plane(self.bary_in(&places[x as usize], far)?),
            plane(self.bary_in(&places[y as usize], far)?),
        );
        for (k, p) in pieces.iter().enumerate() {
            if !alive[k] || dead.contains(&k) || p.face != far {
                continue;
            }
            if [p.a, p.b].contains(&x) || [p.a, p.b].contains(&y) {
                continue;
            }
            let (Some(u), Some(v)) = (
                self.bary_in(&places[p.a as usize], far),
                self.bary_in(&places[p.b as usize], far),
            ) else {
                continue;
            };
            let (u, v) = (plane(u), plane(v));
            let (d1, d2) = (cross2(qx, qy, u), cross2(qx, qy, v));
            let (d3, d4) = (cross2(u, v, qx), cross2(u, v, qy));
            if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                return None;
            }
        }
        let (left, right) = e.regions;
        Some(Shortcut {
            dead,
            piece: Piece {
                a: x,
                b: y,
                face: far,
                left,
                right,
                slack,
            },
            longest,
        })
    }

    /// A place of the seam-cut mesh in space.
    fn position(&self, place: &Place) -> [f64; 3] {
        match *place {
            Place::Face { face, bary } => self.point3(face, bary),
            Place::Edge { edge, t } => {
                let [a, b] = self.flow_mesh.edge_vertices(edge as usize);
                let x = &self.flow_positions;
                [0, 1, 2].map(|k| (1.0 - t) * x[a][k] + t * x[b][k])
            }
        }
    }

    /// The barycentric coordinates of a place in a face of the seam-cut mesh, if it lies on
    /// the closed face.
    fn bary_in(&self, place: &Place, face: u32) -> Option<[f64; 3]> {
        let f = self.flow_mesh.faces()[face as usize];
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
    }

    /// The point of a face with the barycentric coordinates `b`, in space.
    fn point3(&self, face: u32, b: [f64; 3]) -> [f64; 3] {
        let f = self.flow_mesh.faces()[face as usize];
        let x = &self.flow_positions;
        [0, 1, 2].map(|d| (0..3).map(|k| b[k] * x[f[k] as usize][d]).sum())
    }

    /// The longest edge of a face, in space.
    fn longest_edge(&self, face: u32) -> f64 {
        let f = self.flow_mesh.faces()[face as usize];
        let x = &self.flow_positions;
        (0..3)
            .map(|k| {
                let (u, v) = (x[f[k] as usize], x[f[(k + 1) % 3] as usize]);
                (0..3).map(|c| (u[c] - v[c]).powi(2)).sum::<f64>().sqrt()
            })
            .fold(0.0, f64::max)
    }
}

/// Twice the signed area of the triangle `a b c` in the plane.
fn cross2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Whether the rule of `min_offset` moved a point (`place` against where it was `found`).
fn moved_by_the_rule(place: &Place, found: &Place) -> bool {
    match (place, found) {
        (Place::Edge { t, .. }, Place::Edge { t: u, .. }) => t != u,
        (Place::Face { bary, .. }, Place::Face { bary: c, .. }) => bary != c,
        _ => true,
    }
}

/// Whether `q` lies in the convex hull of `points` (closed, to `1e-12`): in one of the
/// triangles of three of them that are not degenerate. With no such triangle every point
/// counts as inside.
fn in_hull(points: &[[f64; 2]], q: [f64; 2]) -> bool {
    let n = points.len();
    let mut any = false;
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let t = [points[i], points[j], points[k]];
                if cross2(t[0], t[1], t[2]).abs() <= 1e-14 {
                    continue;
                }
                any = true;
                if in_closed_triangle(t, q) {
                    return true;
                }
            }
        }
    }
    !any
}

/// Whether `q` lies in the closed triangle `t` (either orientation, to `1e-12`).
fn in_closed_triangle(t: [[f64; 2]; 3], q: [f64; 2]) -> bool {
    let sign = cross2(t[0], t[1], t[2]).signum();
    let eps = 1e-12;
    cross2(t[0], t[1], q) * sign >= -eps
        && cross2(t[1], t[2], q) * sign >= -eps
        && cross2(t[2], t[0], q) * sign >= -eps
}

/// The distance of `p` from the segment `a b` in space.
fn distance_to_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let len2: f64 = ab.iter().map(|v| v * v).sum();
    let s = if len2 > 0.0 {
        (ap.iter().zip(&ab).map(|(u, v)| u * v).sum::<f64>() / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (0..3)
        .map(|k| (ap[k] - s * ab[k]).powi(2))
        .sum::<f64>()
        .sqrt()
}
