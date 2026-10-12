// SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
// SPDX-FileCopyrightText: 2023 Blender Authors
//
// SPDX-License-Identifier: GPL-2.0-or-later

//! Rust port of Blender's face tessellation, so that polygons are split exactly as Blender
//! splits them for `Mesh.loop_triangles`.
//!
//! Ported from Blender v5.0.1 (<https://projects.blender.org/blender/blender>), which is
//! licensed under the GNU General Public License, version 2 or later. **This file is a
//! derivative work of that code and carries the same license**, unlike the rest of the crate.
//! Sources:
//!
//! * `source/blender/blenkernel/intern/mesh_tessellate.cc` — `mesh_calc_tessellation_for_face_impl`
//!   (the per-face rules: triangle, quad with flip test, n-gon via 2D polyfill);
//! * `source/blender/blenlib/intern/polyfill_2d.cc` — the ear clipping triangulator, which is
//!   itself "based on LIBGDX 2013-11-28, APACHE 2.0 licensed" (`EarClippingTriangulator`), with
//!   Blender's changes (even clipping, back-and-forth sweep, convex skip);
//! * `source/blender/blenlib/intern/math_geom.cc` — `is_quad_flip_v3_first_third_fast`,
//!   `axis_dominant_v3_to_m3_negate`, `cross_poly_v2`;
//! * `source/blender/blenlib/intern/math_vector.cc` and `math_vector_inline.cc` —
//!   `ortho_basis_v3v3_v3`, `add_newell_cross_v3_v3v3`, `normalize_v3`;
//! * `source/blender/blenlib/intern/math_matrix_c.cc` — `mul_v2_m3v3`, `transpose_m3`.
//!
//! # Fidelity
//!
//! All arithmetic is done in `f32` with Blender's order of operations, because the result
//! depends on the signs of single precision expressions. Rust never fuses a multiplication and
//! an addition into one FMA instruction; a Blender build whose compiler does (a possibility on
//! ARM64) can therefore decide exact ties differently.
//!
//! One deliberate difference: Blender accelerates the "is any reflex vertex inside this ear"
//! query with a small k-d tree; this port scans the reflex vertices directly. The query result,
//! and hence the triangulation, is identical; only the running time for very large concave
//! n-gons differs (linear instead of logarithmic per query).

type V2 = [f32; 2];
type V3 = [f32; 3];

#[inline]
fn sub_v3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross_v3(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot_v3(a: V3, b: V3) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `is_quad_flip_v3_first_third_fast`: whether the quad `v1 v2 v3 v4` must be split along
/// `v2 v4` instead of the default `v1 v3`, because the two triangles on `v1 v3` would fold
/// over each other.
fn is_quad_flip_v3_first_third_fast(v1: V3, v2: V3, v3: V3, v4: V3) -> bool {
    let d_12 = sub_v3(v2, v1);
    let d_13 = sub_v3(v3, v1);
    let d_14 = sub_v3(v4, v1);
    let cross_a = cross_v3(d_12, d_13);
    let cross_b = cross_v3(d_14, d_13);
    dot_v3(cross_a, cross_b) > 0.0
}

/// The two rows of the matrix that projects a face to 2D: the face's Newell normal, negated,
/// completed to an orthonormal basis (`axis_dominant_v3_to_m3_negate`), of which only the two
/// in-plane axes are needed.
///
/// `corner_positions` holds the positions of the face corners in order.
fn projection_axes(corner_positions: &[V3]) -> [V3; 2] {
    // Newell's method. "Calc normal, flipped: to get a positive 2D cross product."
    let mut normal: V3 = [0.0; 3];
    let mut prev = corner_positions[corner_positions.len() - 1];
    for &curr in corner_positions {
        // add_newell_cross_v3_v3v3
        normal[0] += (prev[1] - curr[1]) * (prev[2] + curr[2]);
        normal[1] += (prev[2] - curr[2]) * (prev[0] + curr[0]);
        normal[2] += (prev[0] - curr[0]) * (prev[1] + curr[1]);
        prev = curr;
    }
    // normalize_v3; a zero-length result falls back to +Z.
    let d = dot_v3(normal, normal);
    if d > 1.0e-35 {
        let d = d.sqrt();
        let f = 1.0 / d;
        normal = [normal[0] * f, normal[1] * f, normal[2] * f];
    } else {
        normal = [0.0, 0.0, 1.0];
    }

    // axis_dominant_v3_to_m3_negate: r_mat[2] = -normal, then ortho_basis_v3v3_v3.
    let n: V3 = [-normal[0], -normal[1], -normal[2]];
    let eps = f32::EPSILON;
    let f = n[0] * n[0] + n[1] * n[1];
    if f > eps {
        let d = 1.0 / f.sqrt();
        let n1: V3 = [n[1] * d, -n[0] * d, 0.0];
        let n2: V3 = [-n[2] * n1[1], n[2] * n1[0], n[0] * n1[1] - n[1] * n1[0]];
        [n1, n2]
    } else {
        // Degenerate case: the normal is (anti)parallel to Z.
        [[if n[2] < 0.0 { -1.0 } else { 1.0 }, 0.0, 0.0], [0.0, 1.0, 0.0]]
    }
}

/// `mul_v2_m3v3` with the transposed axis matrix: the 2D coordinates of `p`.
#[inline]
fn project(axes: &[V3; 2], p: V3) -> V2 {
    [
        axes[0][0] * p[0] + axes[0][1] * p[1] + axes[0][2] * p[2],
        axes[1][0] * p[0] + axes[1][1] * p[1] + axes[1][2] * p[2],
    ]
}

const CONCAVE: i8 = -1;
const TANGENTIAL: i8 = 0;
const CONVEX: i8 = 1;

#[inline]
fn signum_enum(a: f32) -> i8 {
    if a > 0.0 {
        CONVEX
    } else if a == 0.0 {
        TANGENTIAL
    } else {
        CONCAVE
    }
}

/// `area_tri_signed_v2_alt_2x`: twice the signed area of the triangle `v1 v2 v3`.
#[inline]
fn area_tri_signed_v2_alt_2x(v1: V2, v2: V2, v3: V2) -> f32 {
    let d2 = [v2[0] - v1[0], v2[1] - v1[1]];
    let d3 = [v3[0] - v1[0], v3[1] - v1[1]];
    (d2[0] * d3[1]) - (d3[0] * d2[1])
}

#[inline]
fn span_tri_v2_sign(v1: V2, v2: V2, v3: V2) -> i8 {
    signum_enum(area_tri_signed_v2_alt_2x(v3, v2, v1))
}

/// Circular doubly linked list node, stored in an array aligned with the coordinates.
#[derive(Clone, Copy)]
struct PolyIndex {
    next: u32,
    prev: u32,
    index: u32,
    sign: i8,
}

/// Reusable working memory of [`tessellate_face`].
#[derive(Default)]
pub(super) struct PolyFillScratch {
    projected: Vec<V2>,
    indices: Vec<PolyIndex>,
    /// The set Blender keeps in its k-d tree: vertices that were not convex when the fill
    /// started and have neither been clipped nor become convex since.
    reflex: Vec<u32>,
    /// Position of a vertex in `reflex`, or `u32::MAX`.
    reflex_slot: Vec<u32>,
}

struct PolyFill<'a> {
    /// Head of the list (`pf->indices`).
    head: u32,
    indices: &'a mut Vec<PolyIndex>,
    coords: &'a [V2],
    coords_num: u32,
    coords_num_concave: u32,
    reflex: &'a mut Vec<u32>,
    reflex_slot: &'a mut Vec<u32>,
    tris: &'a mut Vec<[u32; 3]>,
}

impl PolyFill<'_> {
    /// `kdtree2d_node_remove`.
    fn reflex_remove(&mut self, index: u32) {
        let slot = self.reflex_slot[index as usize];
        if slot == u32::MAX {
            return;
        }
        self.reflex_slot[index as usize] = u32::MAX;
        let last = self.reflex.pop().expect("slot map and list agree");
        if last != index {
            self.reflex[slot as usize] = last;
            self.reflex_slot[last as usize] = slot;
        }
    }

    /// `kdtree2d_isect_tri`: is any remaining reflex vertex, other than the triangle's own
    /// corners, inside or on the triangle `ind`?
    fn reflex_isect_tri(&self, ind: [u32; 3]) -> bool {
        let vs = ind.map(|i| self.coords[i as usize]);
        let mut min = [f32::MAX; 2];
        let mut max = [-f32::MAX; 2];
        for v in vs {
            for axis in 0..2 {
                if min[axis] > v[axis] {
                    min[axis] = v[axis];
                }
                if max[axis] < v[axis] {
                    max[axis] = v[axis];
                }
            }
        }
        self.reflex.iter().any(|&index| {
            let co = self.coords[index as usize];
            // Bounding box test first, then triangle intersect.
            co[0] >= min[0]
                && co[0] <= max[0]
                && co[1] >= min[1]
                && co[1] <= max[1]
                && span_tri_v2_sign(vs[0], vs[1], co) != CONCAVE
                && span_tri_v2_sign(vs[1], vs[2], co) != CONCAVE
                && span_tri_v2_sign(vs[2], vs[0], co) != CONCAVE
                && index != ind[0]
                && index != ind[1]
                && index != ind[2]
        })
    }

    /// `pf_coord_sign_calc`.
    fn coord_sign_calc(&mut self, pi: u32) {
        let node = self.indices[pi as usize];
        let prev = self.indices[node.prev as usize].index;
        let next = self.indices[node.next as usize].index;
        self.indices[pi as usize].sign = span_tri_v2_sign(
            self.coords[prev as usize],
            self.coords[node.index as usize],
            self.coords[next as usize],
        );
    }

    /// `pf_ear_tip_check`.
    fn ear_tip_check(&self, pi_ear_tip: u32, sign_accept: i8) -> bool {
        // Fast path for circles.
        if self.coords_num_concave == 0 {
            return true;
        }
        let tip = self.indices[pi_ear_tip as usize];
        if tip.sign != sign_accept {
            return false;
        }
        let ind = [tip.index, self.indices[tip.next as usize].index, self.indices[tip.prev as usize].index];
        !self.reflex_isect_tri(ind)
    }

    /// `pf_ear_tip_find`.
    fn ear_tip_find(&self, pi_ear_init: u32, reverse: bool) -> u32 {
        let coords_num = self.coords_num;
        // Two passes: first only convex ears, then accept zero-area (tangential) ones.
        for sign_accept in [CONVEX, TANGENTIAL] {
            let mut pi_ear = pi_ear_init;
            for _ in 0..coords_num {
                if self.ear_tip_check(pi_ear, sign_accept) {
                    return pi_ear;
                }
                let node = self.indices[pi_ear as usize];
                pi_ear = if reverse { node.prev } else { node.next };
            }
        }
        // Desperate mode: no vertex is an ear tip, the polygon is degenerate. Return a convex
        // or tangential vertex if one exists.
        let mut pi_ear = pi_ear_init;
        for _ in 0..coords_num {
            if self.indices[pi_ear as usize].sign != CONCAVE {
                return pi_ear;
            }
            pi_ear = self.indices[pi_ear as usize].next;
        }
        // If all vertices are concave, just return the last one.
        pi_ear
    }

    /// `pf_ear_tip_cut` and `pf_coord_remove`.
    fn ear_tip_cut(&mut self, pi_ear_tip: u32) {
        let tip = self.indices[pi_ear_tip as usize];
        self.tris.push([self.indices[tip.prev as usize].index, tip.index, self.indices[tip.next as usize].index]);

        if !self.reflex.is_empty() {
            self.reflex_remove(tip.index);
        }
        self.indices[tip.next as usize].prev = tip.prev;
        self.indices[tip.prev as usize].next = tip.next;
        if self.head == pi_ear_tip {
            self.head = tip.next;
        }
        self.coords_num -= 1;
    }

    /// `pf_triangulate`.
    fn triangulate(&mut self) {
        let mut pi_ear_init = self.head;
        let mut reverse = false;

        while self.coords_num > 3 {
            let pi_ear = self.ear_tip_find(pi_ear_init, reverse);
            if self.indices[pi_ear as usize].sign != CONVEX {
                self.coords_num_concave -= 1;
            }
            let pi_prev = self.indices[pi_ear as usize].prev;
            let pi_next = self.indices[pi_ear as usize].next;

            self.ear_tip_cut(pi_ear);

            // The type of the two vertices adjacent to the clipped vertex may have changed.
            for pi in [pi_prev, pi_next] {
                if self.indices[pi as usize].sign != CONVEX {
                    self.coord_sign_calc(pi);
                    if self.indices[pi as usize].sign == CONVEX {
                        self.coords_num_concave -= 1;
                        self.reflex_remove(self.indices[pi as usize].index);
                    }
                }
            }

            pi_ear_init = if reverse {
                self.indices[pi_prev as usize].prev
            } else {
                self.indices[pi_next as usize].next
            };
            if self.indices[pi_ear_init as usize].sign != CONVEX {
                // Take the extra step since this ear isn't a good candidate.
                let node = self.indices[pi_ear_init as usize];
                pi_ear_init = if reverse { node.prev } else { node.next };
                reverse = !reverse;
            }
        }

        if self.coords_num == 3 {
            let a = self.indices[self.head as usize];
            let b = self.indices[a.next as usize];
            let c = self.indices[b.next as usize];
            self.tris.push([a.index, b.index, c.index]);
        }
    }
}

/// `BLI_polyfill_calc` with `coords_sign = 1`: triangulates the polygon `coords` (which the
/// caller projected so that it winds the way Blender expects) and appends `coords.len() - 2`
/// triangles of indices into `coords` to `tris`.
///
/// `coords.len()` must be at least 3.
fn polyfill_calc(
    coords: &[V2],
    indices: &mut Vec<PolyIndex>,
    reflex: &mut Vec<u32>,
    reflex_slot: &mut Vec<u32>,
    tris: &mut Vec<[u32; 3]>,
) {
    let coords_num = coords.len() as u32;

    // polyfill_prepare
    indices.clear();
    indices.extend((0..coords_num).map(|i| PolyIndex {
        next: if i + 1 == coords_num { 0 } else { i + 1 },
        prev: if i == 0 { coords_num - 1 } else { i - 1 },
        index: i,
        sign: CONVEX,
    }));
    reflex.clear();
    reflex_slot.clear();
    reflex_slot.resize(coords_num as usize, u32::MAX);

    let mut pf = PolyFill {
        head: 0,
        indices,
        coords,
        coords_num,
        coords_num_concave: 0,
        reflex,
        reflex_slot,
        tris,
    };
    for i in 0..coords_num {
        pf.coord_sign_calc(i);
        if pf.indices[i as usize].sign != CONVEX {
            pf.coords_num_concave += 1;
        }
    }

    // polyfill_calc: kdtree2d_init collects every vertex that is not convex.
    if pf.coords_num_concave != 0 {
        for i in 0..coords_num {
            if pf.indices[i as usize].sign != CONVEX {
                pf.reflex_slot[i as usize] = pf.reflex.len() as u32;
                pf.reflex.push(i);
            }
        }
    }

    pf.triangulate();
}

/// `mesh_calc_tessellation_for_face_impl` (without a precomputed normal): splits one face,
/// given the positions of its corners in order, into `corner_positions.len() - 2` triangles.
/// The triangles are appended to `tris` as indices into `corner_positions`, in Blender's order
/// and with Blender's corner order.
///
/// `corner_positions.len()` must be at least 3.
pub(super) fn tessellate_face(corner_positions: &[V3], scratch: &mut PolyFillScratch, tris: &mut Vec<[u32; 3]>) {
    match corner_positions.len() {
        3 => tris.push([0, 1, 2]),
        4 => {
            let p = corner_positions;
            if is_quad_flip_v3_first_third_fast(p[0], p[1], p[2], p[3]) {
                // Flip out of degenerate 0-2 state.
                tris.push([0, 1, 3]);
                tris.push([1, 2, 3]);
            } else {
                tris.push([0, 1, 2]);
                tris.push([0, 2, 3]);
            }
        }
        _ => {
            let axes = projection_axes(corner_positions);
            let PolyFillScratch { projected, indices, reflex, reflex_slot } = scratch;
            projected.clear();
            projected.extend(corner_positions.iter().map(|&p| project(&axes, p)));
            polyfill_calc(projected, indices, reflex, reflex_slot, tris);
        }
    }
}
