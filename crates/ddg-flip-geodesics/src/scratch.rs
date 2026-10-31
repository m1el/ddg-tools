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

    /// The angle of a wedge (Eq. 1, [`crate::wedge_angle`]): its corner angles summed in the
    /// order given, bitwise as `wedge_angle`.
    #[inline]
    pub(crate) fn wedge_angle<T: Triangulation>(&mut self, tri: &T, corners: &[u32]) -> f64 {
        corners.iter().map(|&h| self.angle(tri, h as usize)).sum()
    }
}

/// A joint found in the star of its vertex: the positions of its two segment ends in
/// [`VertexStar::entries`] (`p` the end of `ab`, `q` the start of `bc`; SPEC §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Located {
    pub(crate) joint: Joint,
    pub(crate) vertex: usize,
    pub(crate) p: usize,
    pub(crate) q: usize,
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
        self.star.rebuild(complex, network, v);
        let p = self.star.position(joint.incoming, false).ok_or_else(bad)?;
        let q = self.star.position(joint.outgoing, true).ok_or_else(bad)?;
        Ok(Located {
            joint,
            vertex: v,
            p,
            q,
        })
    }

    /// The sweep of the wedge on `side` into [`corners`](Self::corners), in sweep order (the
    /// right wedge from `p` to `q`, the left one from `q` to `p`); whether it crosses the
    /// boundary gap. `at` must be the last joint located.
    fn sweep(&mut self, at: Located, side: Side) -> bool {
        let (from, to) = match side {
            Side::Right => (at.p, at.q),
            Side::Left => (at.q, at.p),
        };
        self.star.sweep_corners(from, to, &mut self.corners)
    }

    /// The wedge angle on `side` (bitwise [`JointGeometry::angle`](crate::JointGeometry::angle)).
    /// `at` must be the last joint located.
    pub(crate) fn angle<T: Triangulation>(&mut self, tri: &T, at: Located, side: Side) -> f64 {
        if self.sweep(at, side) {
            f64::INFINITY
        } else {
            self.angles.wedge_angle(tri, &self.corners)
        }
    }

    /// The corners of the wedge on `side` ordered from `ba` to `bc` into
    /// [`corners`](Self::corners)
    /// ([`JointGeometry::corners_from_incoming`](crate::JointGeometry::corners_from_incoming)).
    /// `at` must be the last joint located.
    pub(crate) fn corners_from_incoming(&mut self, at: Located, side: Side) {
        self.sweep(at, side);
        if side == Side::Left {
            self.corners.reverse();
        }
    }

    /// Whether the wedge on `side` is flexible: no segment end strictly inside it
    /// ([`JointGeometry::is_flexible`](crate::JointGeometry::is_flexible)). `at` must be the
    /// last joint located.
    pub(crate) fn is_flexible(&self, at: Located, side: Side) -> bool {
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
