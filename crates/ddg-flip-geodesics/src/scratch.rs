//! Reusable working memory of the shortening (SPEC §4.1, §4.2): the star of the joint's
//! vertex, the corners of a wedge, the buffers of a replacement, and a memo of corner angles.
//! Nothing here changes a result: every value is the one [`JointGeometry`] and
//! [`Triangulation::corner_angle`] give, computed in the same order (F2), without allocating
//! per call.
//!
//! [`JointGeometry`]: crate::JointGeometry

use std::f64::consts::PI;

use ddg_intrinsic::Triangulation;

use crate::error::FlipError;
use crate::network::{Joint, Network, Replaced, Segment};
use crate::wedge::{Side, VertexStar};

/// A memo of corner angles (F2's cache). One slot per face, direct mapped by a hash of the face
/// id; a slot holds the bits of the face's three lengths (by slot, as
/// [`Triangulation::face_lengths`]) and the angles of its corners computed so far with
/// [`Triangulation::corner_angle`]. A slot is used only when the face's lengths are bitwise the
/// ones it was filled with; a corner angle is a function of those three lengths, so the memo
/// returns exactly what `corner_angle` returns now, whatever flipped since. Colliding faces
/// evict each other (only time is lost).
#[derive(Debug, Clone)]
pub(crate) struct AngleCache {
    slots: Vec<Slot>,
    shift: u32,
}

#[derive(Debug, Clone, Copy, Default)]
struct Slot {
    /// The face id plus one; `0` for an empty slot.
    face: u32,
    lengths: [u64; 3],
    /// `0.0` for a corner not computed yet (a corner angle is `> 0`).
    angles: [f64; 3],
}

impl AngleCache {
    /// A memo with `2^bits` slots.
    pub(crate) fn new(bits: u32) -> AngleCache {
        AngleCache {
            slots: vec![Slot::default(); 1 << bits],
            shift: 64 - bits,
        }
    }

    /// [`Triangulation::corner_angle`] of corner `h`, bitwise.
    #[inline]
    pub(crate) fn angle<T: Triangulation>(&mut self, tri: &T, h: usize) -> f64 {
        let f = h / 3;
        let lengths = tri.face_lengths(f).map(f64::to_bits);
        let key = f as u32 + 1;
        let i = (u64::from(key).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize;
        let slot = &mut self.slots[i];
        if slot.face != key || slot.lengths != lengths {
            *slot = Slot {
                face: key,
                lengths,
                angles: [0.0; 3],
            };
        }
        let a = &mut slot.angles[h % 3];
        if *a == 0.0 {
            *a = tri.corner_angle(h);
        }
        *a
    }

    /// The memoized angle of corner `h` of face `f`, whose lengths have the bits `lengths`.
    #[inline]
    fn peek(&self, f: usize, lengths: [u64; 3], k: usize) -> Option<f64> {
        let key = f as u32 + 1;
        let i = (u64::from(key).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize;
        let slot = &self.slots[i];
        (slot.face == key && slot.lengths == lengths && slot.angles[k] != 0.0)
            .then_some(slot.angles[k])
    }

    /// Whether `θ_g + θ_h < π` for the corner angles of corners `g` and `h`
    /// ([`Triangulation::corner_angle`], summed in floating point): bitwise the comparison
    /// `angle(g) + angle(h) < PI`, mostly without computing the angles.
    ///
    /// With Kahan's factors of a corner (`ddg_core::geometry::corner_angle`:
    /// `θ = 2 atan2(√num, √den)`), `T = tan(θ/2) = √(num/den)` and, for two corners with
    /// `P = T_g T_h`, `θ_g + θ_h − π = 2 atan((P − 1)/(T_g + T_h))`. So when
    /// `|P² − 1| > (P + 1)(T_g + T_h) δ` with `δ = 1e-13` the sign of `P² − 1` decides
    /// `θ_g + θ_h < π` exactly, and also the floating-point sum (its error is below `1e-14`:
    /// a few ulps of each angle and of the sum, and `PI` is within `1.3e-16` of `π`). The bound
    /// uses `T ≤ max(1, T²)`; the products carry a relative error below `1e-15`, covered by an
    /// extra `1e-14 P²`. When both angles are memoized, or the test is inconclusive, the angles
    /// are summed (and memoized).
    pub(crate) fn sum_below_pi<T: Triangulation>(&mut self, tri: &T, g: usize, h: usize) -> bool {
        let (lg, lh) = (tri.face_lengths(g / 3), tri.face_lengths(h / 3));
        if let (Some(a), Some(b)) = (
            self.peek(g / 3, lg.map(f64::to_bits), g % 3),
            self.peek(h / 3, lh.map(f64::to_bits), h % 3),
        ) {
            return a + b < PI;
        }
        let (ng, dg) = kahan_factors(lg, g % 3);
        let (nh, dh) = kahan_factors(lh, h % 3);
        let (sg, sh) = (ng / dg, nh / dh);
        let q = sg * sh;
        let bound = (q.max(1.0) + 1.0) * (sg.max(1.0) + sh.max(1.0)) * 1e-13 + 1e-14 * q;
        let decided = if q < 1.0 - bound {
            Some(true)
        } else if q > 1.0 + bound {
            Some(false)
        } else {
            None
        };
        match decided {
            Some(below) => {
                debug_assert_eq!(
                    below,
                    self.angle(tri, g) + self.angle(tri, h) < PI,
                    "the filter of beta < pi disagrees with the angles"
                );
                below
            }
            None => self.angle(tri, g) + self.angle(tri, h) < PI,
        }
    }

    /// The angle of a wedge (Eq. 1, [`crate::wedge_angle`]): its corner angles summed in the
    /// order given, bitwise as `wedge_angle`.
    #[inline]
    pub(crate) fn wedge_angle<T: Triangulation>(&mut self, tri: &T, corners: &[u32]) -> f64 {
        corners.iter().map(|&h| self.angle(tri, h as usize)).sum()
    }
}

/// Kahan's factors `(num, den)` of the corner at slot `k` of a face with lengths `l` (by
/// slot), `θ = 2 atan2(√num, √den)`, as `ddg_core::geometry::corner_angle` computes them for
/// [`Triangulation::corner_angle`] (opposite side `l[k + 1]`, sides `l[k]`, `l[k + 2]`). Only
/// used to decide comparisons ([`AngleCache::sum_below_pi`]), whose margin does not need them
/// bitwise.
fn kahan_factors(l: [f64; 3], k: usize) -> (f64, f64) {
    let (c, s1, s2) = (l[(k + 1) % 3], l[k], l[(k + 2) % 3]);
    let (a, b) = if s1 >= s2 { (s1, s2) } else { (s2, s1) };
    let mu = if b >= c { c - (a - b) } else { b - (a - c) };
    (((a - b) + c) * mu, (a + (b + c)) * ((a - c) + b))
}

/// A joint found at its vertex `b`. Either through the star of `b`: the positions `p` (the
/// end of `ab`) and `q` (the start of `bc`) in [`VertexStar::entries`] (SPEC §3); or, when `b`
/// is interior and the two segment ends lie on different edge ends (`fast`), by the halfedges
/// `hp`, `hq` leaving `b` along `ab` and `bc`: the wedges are then walked around `b` directly
/// (`h → twin(prev(h))`, counter-clockwise), which meets the same corners in the same order as
/// the star's sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Located {
    pub(crate) joint: Joint,
    pub(crate) vertex: usize,
    pub(crate) p: usize,
    pub(crate) q: usize,
    pub(crate) fast: bool,
    pub(crate) hp: u32,
    pub(crate) hq: u32,
}

/// The working memory of one shortening run.
#[derive(Debug, Clone)]
pub(crate) struct Scratch {
    /// The star of the vertex last located.
    pub(crate) star: VertexStar,
    /// The corners of the wedge last swept.
    pub(crate) corners: Vec<u32>,
    /// Edges of a wedge, for the Δ-complex cases (App. A).
    pub(crate) edges: Vec<usize>,
    /// The outer arc of the last Alg. 1 call.
    pub(crate) new: Vec<(Segment, usize)>,
    /// What the last replacement did.
    pub(crate) replaced: Replaced,
    pub(crate) angles: AngleCache,
}

impl Scratch {
    /// Working memory with a memo of `2^bits` faces.
    pub(crate) fn new(bits: u32) -> Scratch {
        Scratch {
            star: VertexStar {
                vertex: 0,
                boundary: false,
                ends: Vec::new(),
                entries: Vec::new(),
                first_entry: Vec::new(),
            },
            corners: Vec::new(),
            edges: Vec::new(),
            new: Vec::new(),
            replaced: Replaced::default(),
            angles: AngleCache::new(bits),
        }
    }

    /// Builds the star of the joint's vertex and finds the joint in it, as
    /// [`JointGeometry::new`](crate::JointGeometry::new) does.
    ///
    /// # Errors
    /// [`FlipError::IndexOutOfRange`] if the joint is not two consecutive live segments.
    pub(crate) fn locate<T: Triangulation>(
        &mut self,
        tri: &T,
        network: &Network,
        joint: Joint,
    ) -> Result<Located, FlipError> {
        let bad = || FlipError::IndexOutOfRange {
            what: "joint (live consecutive segments)",
            index: joint.outgoing.0 as usize,
            bound: network.segment_id_bound(),
        };
        let outgoing = network.segment(joint.outgoing).ok_or_else(bad)?;
        if outgoing.prev != Some(joint.incoming) {
            return Err(bad());
        }
        let complex = tri.complex();
        let v = outgoing.segment.tail(complex);
        if !complex.is_boundary_vertex(v)
            && let Some(incoming) = network.segment(joint.incoming)
            && incoming.segment.head(complex) == v
            && let (Some(hp), Some(hq)) = (
                incoming.segment.reversed().halfedge(complex),
                outgoing.segment.halfedge(complex),
            )
            && hp != hq
        {
            let at = Located {
                joint,
                vertex: v,
                p: 0,
                q: 0,
                fast: true,
                hp: hp as u32,
                hq: hq as u32,
            };
            #[cfg(debug_assertions)]
            self.check_fast(tri, network, at)?;
            return Ok(at);
        }
        self.locate_star(tri, network, joint)
    }

    /// [`locate`](Self::locate) through the star of the joint's vertex.
    ///
    /// # Errors
    /// As [`locate`](Self::locate).
    pub(crate) fn locate_star<T: Triangulation>(
        &mut self,
        tri: &T,
        network: &Network,
        joint: Joint,
    ) -> Result<Located, FlipError> {
        let bad = || FlipError::IndexOutOfRange {
            what: "joint (live consecutive segments)",
            index: joint.outgoing.0 as usize,
            bound: network.segment_id_bound(),
        };
        let outgoing = network.segment(joint.outgoing).ok_or_else(bad)?;
        if outgoing.prev != Some(joint.incoming) {
            return Err(bad());
        }
        let complex = tri.complex();
        let v = outgoing.segment.tail(complex);
        self.star.rebuild(complex, network, v);
        let p = self.star.position(joint.incoming, false).ok_or_else(bad)?;
        let q = self.star.position(joint.outgoing, true).ok_or_else(bad)?;
        Ok(Located {
            joint,
            vertex: v,
            p,
            q,
            fast: false,
            hp: 0,
            hq: 0,
        })
    }

    /// Debug builds: the walk of a fast joint gives bitwise the star's angles, corners and
    /// flexibility.
    #[cfg(debug_assertions)]
    fn check_fast<T: Triangulation>(
        &mut self,
        tri: &T,
        network: &Network,
        at: Located,
    ) -> Result<(), FlipError> {
        let star = self.locate_star(tri, network, at.joint)?;
        for side in [Side::Left, Side::Right] {
            self.corners_from_incoming(tri, star, side);
            let by_star = self.corners.clone();
            let angle_star = self.angle(tri, star, side).to_bits();
            let flexible_star = self.is_flexible(network, tri, star, side);
            self.corners_from_incoming(tri, at, side);
            assert_eq!(
                by_star, self.corners,
                "walked corners differ from the star's"
            );
            assert_eq!(angle_star, self.angle(tri, at, side).to_bits());
            assert_eq!(flexible_star, self.is_flexible(network, tri, at, side));
        }
        Ok(())
    }

    /// The corners walked counter-clockwise around the vertex from halfedge `from` (leaving
    /// it) up to, not including, `to`, into [`corners`](Self::corners).
    fn walk<T: Triangulation>(&mut self, tri: &T, from: u32, to: u32) {
        let complex = tri.complex();
        self.corners.clear();
        let mut h = from as usize;
        loop {
            self.corners.push(h as u32);
            h = complex
                .twin(complex.prev(h))
                .expect("an interior vertex has no boundary gap");
            if h == to as usize {
                break;
            }
        }
    }

    /// The sweep of the wedge on `side` into [`corners`](Self::corners), in sweep order (the
    /// right wedge from `p` to `q`, the left one from `q` to `p`); whether it crosses the
    /// boundary gap. `at` must be the last joint located.
    fn sweep<T: Triangulation>(&mut self, tri: &T, at: Located, side: Side) -> bool {
        if at.fast {
            match side {
                Side::Right => self.walk(tri, at.hp, at.hq),
                Side::Left => self.walk(tri, at.hq, at.hp),
            }
            return false;
        }
        let (from, to) = match side {
            Side::Right => (at.p, at.q),
            Side::Left => (at.q, at.p),
        };
        self.star.sweep_corners(from, to, &mut self.corners)
    }

    /// The wedge angle on `side` (bitwise [`JointGeometry::angle`](crate::JointGeometry::angle)).
    /// `at` must be the last joint located.
    pub(crate) fn angle<T: Triangulation>(&mut self, tri: &T, at: Located, side: Side) -> f64 {
        if self.sweep(tri, at, side) {
            f64::INFINITY
        } else {
            self.angles.wedge_angle(tri, &self.corners)
        }
    }

    /// The corners of the wedge on `side` ordered from `ba` to `bc` into
    /// [`corners`](Self::corners)
    /// ([`JointGeometry::corners_from_incoming`](crate::JointGeometry::corners_from_incoming)).
    /// `at` must be the last joint located.
    pub(crate) fn corners_from_incoming<T: Triangulation>(
        &mut self,
        tri: &T,
        at: Located,
        side: Side,
    ) {
        self.sweep(tri, at, side);
        if side == Side::Left {
            self.corners.reverse();
        }
    }

    /// Whether the wedge on `side` is flexible: no segment end strictly inside it
    /// ([`JointGeometry::is_flexible`](crate::JointGeometry::is_flexible)). `at` must be the
    /// last joint located.
    ///
    /// For a fast joint (as `gc: fg.cpp:589-658`): the joint's segments are outermost on the
    /// wedge side of their edges, and no edge strictly inside the wedge carries a segment.
    /// Within an edge end at `b`, the segment met last turning counter-clockwise is the one
    /// outermost on the side of the face of the halfedge leaving `b`
    /// ([`Network::outermost`]), the one met first is outermost on the side of its twin.
    pub(crate) fn is_flexible<T: Triangulation>(
        &self,
        network: &Network,
        tri: &T,
        at: Located,
        side: Side,
    ) -> bool {
        if at.fast {
            let complex = tri.complex();
            let (from, to, first, last) = match side {
                Side::Right => (at.hp, at.hq, at.joint.incoming, at.joint.outgoing),
                Side::Left => (at.hq, at.hp, at.joint.outgoing, at.joint.incoming),
            };
            let twin = |h: u32| {
                complex
                    .twin(h as usize)
                    .expect("an interior vertex has no boundary gap")
            };
            if network.outermost(complex, from as usize) != Some(first)
                || network.outermost(complex, twin(to)) != Some(last)
            {
                return false;
            }
            let mut h = twin(complex.prev(from as usize) as u32);
            while h != to as usize {
                if network.carries_segment(complex.edge(h)) {
                    return false;
                }
                h = twin(complex.prev(h) as u32);
            }
            return true;
        }
        let m = self.star.entries.len();
        let (from, to) = match side {
            Side::Right => (at.p, at.q),
            Side::Left => (at.q, at.p),
        };
        (to + m - from) % m <= 1
    }

    /// Both wedge angles `(left, right)` of the last joint located.
    pub(crate) fn angles<T: Triangulation>(&mut self, tri: &T, at: Located) -> (f64, f64) {
        let right = self.angle(tri, at, Side::Right);
        let left = self.angle(tri, at, Side::Left);
        (left, right)
    }
}

/// The sides to queue ([`JointGeometry::sides_to_shorten`](crate::JointGeometry::sides_to_shorten)):
/// none if straight; the smaller side, and the larger one too if it is also below
/// `π − tolerance`.
pub(crate) fn sides_to_shorten(left: f64, right: f64, tolerance: f64) -> ([(Side, f64); 2], usize) {
    let (min, max) = if left < right {
        ((Side::Left, left), (Side::Right, right))
    } else {
        ((Side::Right, right), (Side::Left, left))
    };
    // `take_while(α ≤ π − tol)` over (min, max), as `sides_to_shorten`.
    let bent = |a: f64| a <= PI - tolerance;
    let n = if !bent(min.1) {
        0
    } else if !bent(max.1) {
        1
    } else {
        2
    };
    ([min, max], n)
}
