//! Insertion of a vertex in a face (paper Sec. 3.3, Supp. Alg. 5, 6; SPEC §4.4, I6, I15, I25).
//!
//! A vertex `m` is inserted at the barycentric point `v` of face `ijk` of `T¹` and joined to the
//! three corners. Everything is decided before anything changes:
//!
//! 1. **Lengths and refusal** (I25): [`plan_face_split`], the same function the signpost
//!    representation calls, so both refuse the same insertions and give the new edges bitwise
//!    the same lengths (`it-point-distances`; the paper's Eq. 5, `ic-tangent-length`, is the same
//!    quantity, proved equal in Lean, and is not evaluated).
//! 2. **The arcs of the face.** Each side's crossings come from
//!    [`extract_intrinsic_edge`](IntegerTriangulation::extract_intrinsic_edge) (WP2,
//!    `O(n_e + walk)` per side): the input edge of every crossing, exact, and its position
//!    `distance / ℓ` along the side. The arcs inside the face are given by the integers alone
//!    ([`FaceArrangement`]): `c_k` arcs around each corner `k` (the `t`-th from crossing `t` of
//!    side `k` to crossing `n_{k+2} − 1 − t` of side `k + 2`) and, when a corner `κ` emanates,
//!    `e` arcs from `κ` to crossing `c_{κ+1} + s` of the opposite side. Each arc is the straight
//!    chord between its two end points in the face laid out in the plane.
//! 3. **The region** (I15). The regions of the arrangement are enumerated combinatorially in
//!    the oracle's labels ([`Region`]: `corner k t`, `middle s`); every region is the triangle
//!    cut by at most five chords, so each has a *margin*: the smallest signed distance of `m`
//!    to its bounding chords, positive inside. Floating point only **picks** the region: the one
//!    with the largest margin, ties by the oracle's order (the region of `m` in exact
//!    arithmetic is the only one with a positive margin). The new normal coordinates are then
//!    computed from that region exactly ([`FaceArrangement::split_counts`],
//!    `ic-face-split-normal`: the number of arcs that separate the region from each corner; Lean
//!    `faceSplit?_gapRegion`), so the result is valid whichever region is picked (Lean
//!    `splitCount_faceValid`). This replaces the paper's slack rounding of Alg. 5 lines 10–16 and
//!    Alg. 6 lines 11–14, and the reference's snapping (SPEC I15).
//! 4. **I6.** The distance of `m` to the nearest arc (the nearest bounding chord of its region,
//!    as a segment), in units of the face's shortest side, is reported
//!    ([`RegionChoice::curve_distance`]); with [`FaceSplitOptions::too_close`] above it the
//!    insertion is refused with [`IntegerError::TooCloseToCurve`] and nothing changes. The
//!    default `0` never refuses.
//! 5. **The roundabouts** of the new halfedges `a → m` with `a ∈ V⁰` (Eq. 4's analogue, formula
//!    id `ic-roundabout-update`, tier C, not tagged): with `a` the tail of side `k`, the input
//!    halfedges at `a` between side `k` and `a → m` are the one along side `k` (`n⁻`) and the
//!    arcs that `a` emanates on that side of `m` (`s` for the fan region `middle s`, `e` for a
//!    corner region of `κ + 2`, `0` otherwise): `r_{am} = r_{side k} + n⁻ + those`, reduced as
//!    [`LocalIndices::modulus`](crate::LocalIndices::modulus) says.
//! 6. **The place** `q⁰_m`: `RecoverBarycentric` (App. A, [`crate::recover`]) over every corner
//!    of the region (SPEC I16).
//!
//! Then the connectivity changes by [`Complex::split_face`](ddg_intrinsic::Complex::split_face)
//! (see [`ddg_intrinsic::complex::FaceSplit`] for the ids). The point query of SPEC §4.9 is the same
//! without the change ([`IntegerTriangulation::locate_in_face`],
//! [`IntegerTriangulation::point_on_input`]).

use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic::{
    EdgeTrace, EdgeTraceOptions, Insertion, IntrinsicError, SurfacePoint, plan_face_split,
};

use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{FaceCoordinates, n_minus, roundabout_update};
use crate::geometry::StripReport;
use crate::recover::{RecoverReport, RecoveredPoint, recover_in_input_face};

/// A region of the arrangement of the arcs of a face, in the Lean oracle's labels
/// (`lean/oracle/README.md`, tasks `regions`, `split`; Lean `Ddg.IntegerCoords.Region`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Region {
    /// `corner k t` (`t < c_k`): the region between the `(t − 1)`-st and the `t`-th arc around
    /// corner `k`, counted from the corner (`t = 0` contains the corner).
    Corner {
        /// The corner `k` (`0, 1, 2`).
        corner: usize,
        /// `t`: the number of arcs around the corner between the region and the corner.
        index: i64,
    },
    /// `middle s`: what is left: the central region (`s = 0`) when no corner emanates, or the
    /// `s`-th region of the fan (`0 ≤ s ≤ e`, counted from corner `κ + 1`) when corner `κ`
    /// emanates `e > 0` arcs.
    Middle {
        /// `s`.
        index: i64,
    },
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Region::Corner { corner, index } => write!(f, "corner {corner} {index}"),
            Region::Middle { index } => write!(f, "middle {index}"),
        }
    }
}

/// The arcs of a valid face, decoded from its three normal coordinates (Eq. 1, 2): the corner
/// coordinates and the emanating corner `κ` with its count `e` (`κ = 0`, `e = 0` when nothing
/// emanates), as the oracle's `decode?`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceArrangement {
    /// `n⁺`, `c`, `e` of the face.
    pub coords: FaceCoordinates,
    /// The emanating corner `κ` (`0` when none emanates).
    pub kappa: usize,
    /// `e_κ` (`0` when none emanates).
    pub fan: i64,
}

impl FaceArrangement {
    /// The arrangement of the face with side counts `n` (side `k` from corner `k` to `k + 1`),
    /// or `None` for a face that is not valid (an odd numerator of Eq. 1).
    ///
    /// [formula: ic-face-decoding]
    pub fn of(n: [i64; 3]) -> Option<FaceArrangement> {
        let coords = FaceCoordinates::of(n).ok()?;
        let (kappa, fan) = (0..3)
            .find(|&k| coords.e[k] > 0)
            .map_or((0, 0), |k| (k, coords.e[k]));
        Some(FaceArrangement { coords, kappa, fan })
    }

    /// Whether `region` is one of this face's.
    pub fn contains(&self, region: Region) -> bool {
        match region {
            Region::Corner { corner, index } => {
                corner < 3 && (0..self.coords.c[corner]).contains(&index)
            }
            Region::Middle { index } => (0..=self.fan).contains(&index),
        }
    }

    /// The number of regions: `c₀ + c₁ + c₂ + e + 1`.
    pub fn n_regions(&self) -> u128 {
        self.coords.c.iter().map(|&c| c as u128).sum::<u128>() + self.fan as u128 + 1
    }

    /// The regions in the oracle's order: the corner regions of corners `0, 1, 2` by `t`, then
    /// `middle 0 … e`.
    ///
    /// [formula: ic-face-split-normal]
    pub fn regions(&self) -> impl Iterator<Item = Region> + '_ {
        (0..3)
            .flat_map(move |k| {
                (0..self.coords.c[k]).map(move |t| Region::Corner {
                    corner: k,
                    index: t,
                })
            })
            .chain((0..=self.fan).map(|s| Region::Middle { index: s }))
    }

    /// The normal coordinates `n_{m0}, n_{m1}, n_{m2}` of the edges from a vertex `m` inserted
    /// in `region` to the corners `0, 1, 2`: the number of arcs that separate the region from
    /// each corner (an arc that ends at the corner does not). Corner region `(k, t)`: `t` arcs to
    /// corner `k`, and `c_j + c_k − t` to each other corner `j`, plus the `e` fan arcs when
    /// `{k, j} = {κ + 1, κ + 2}`. Middle region `s`: `c_j`, plus `s` (for `j = κ + 1`) or `e − s`
    /// (for `j = κ + 2`). The oracle's `splitCount`; Lean proves it is the separation count of
    /// every region (`faceSplit?_gapRegion`) and that the three new faces are valid
    /// (`splitCount_faceValid`). `None` for a region the face does not have, or a count beyond
    /// `i64`.
    ///
    /// [formula: ic-face-split-normal]
    pub fn split_counts(&self, region: Region) -> Option<[i64; 3]> {
        if !self.contains(region) {
            return None;
        }
        let c = |j: usize| i128::from(self.coords.c[j]);
        let (next, prev) = ((self.kappa + 1) % 3, (self.kappa + 2) % 3);
        let e = i128::from(self.fan);
        let mut out = [0i64; 3];
        for (j, slot) in out.iter_mut().enumerate() {
            let value = match region {
                Region::Corner { corner: k, index } => {
                    let t = i128::from(index);
                    if j == k {
                        t
                    } else {
                        let fan = if (k == next && j == prev) || (k == prev && j == next) {
                            e
                        } else {
                            0
                        };
                        c(j) + c(k) - t + fan
                    }
                }
                Region::Middle { index } => {
                    let s = i128::from(index);
                    if j == next {
                        c(j) + s
                    } else if j == prev {
                        c(j) + e - s
                    } else {
                        c(j)
                    }
                }
            };
            *slot = i64::try_from(value).ok()?;
        }
        Some(out)
    }

    /// The number of arcs that corner `k` emanates between its side `k` (towards corner
    /// `k + 1`) and a vertex inserted in `region`, counter-clockwise around the corner: the
    /// input halfedges leaving the corner inside the face before the new edge.
    pub fn emanating_before(&self, k: usize, region: Region) -> i64 {
        if self.fan == 0 || k != self.kappa {
            return 0;
        }
        match region {
            Region::Middle { index } => index,
            Region::Corner { corner, .. } if corner == (self.kappa + 2) % 3 => self.fan,
            Region::Corner { .. } => 0,
        }
    }
}

/// The regions of the face with side counts `n` in the oracle's order (task `regions`), or
/// `None` for a face that is not valid.
///
/// [formula: ic-face-split-normal]
pub fn face_regions(n: [i64; 3]) -> Option<Vec<Region>> {
    FaceArrangement::of(n).map(|a| a.regions().collect())
}

/// The face-split counts of `region` of the face with side counts `n` (task `split`):
/// [`FaceArrangement::split_counts`]; `None` for a face that is not valid or a region it does
/// not have.
///
/// [formula: ic-face-split-normal]
pub fn split_counts(n: [i64; 3], region: Region) -> Option<[i64; 3]> {
    FaceArrangement::of(n)?.split_counts(region)
}

/// Options of [`IntegerTriangulation::split_face`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceSplitOptions {
    /// SPEC I6: refuse the insertion ([`IntegerError::TooCloseToCurve`], nothing changed) when
    /// the point is closer than this to an arc of the face (an input edge crossing it), in
    /// units of the face's shortest side. `0` (the default) never refuses. Callers that need the
    /// inserted point exactly where they asked (refinement, SPEC §5) set it.
    pub too_close: f64,
}

impl Default for FaceSplitOptions {
    fn default() -> Self {
        FaceSplitOptions { too_close: 0.0 }
    }
}

/// How the region of a point was chosen (SPEC I6, I15).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegionChoice {
    /// The region.
    pub region: Region,
    /// Its margin: the smallest signed distance of the point to the chords that bound the
    /// region, positive inside, in units of the face's shortest side; `+∞` for a face without
    /// arcs. Negative only when rounding put the point outside every region (then the region
    /// with the largest margin is the one it is nearest to).
    pub margin: f64,
    /// The region with the second largest margin and that margin, if the face has more than
    /// one region. In exact arithmetic it is negative.
    pub runner_up: Option<(Region, f64)>,
    /// SPEC I6: the distance of the point to the nearest arc of the face (the nearest bounding
    /// chord of its region, as a segment), in units of the face's shortest side; `None` for a
    /// face without arcs.
    pub curve_distance: Option<f64>,
    /// The number of regions of the face.
    pub n_regions: u128,
}

/// Where a point of a face of the triangulation lies with respect to the input mesh: its
/// region, the counts a vertex inserted there gets, and its place recovered by App. A
/// ([`IntegerTriangulation::locate_in_face`]).
#[derive(Debug, Clone, PartialEq)]
pub struct FaceLocation {
    /// The face of the triangulation.
    pub face: u32,
    /// The arrangement of its arcs.
    pub arrangement: FaceArrangement,
    /// The region of the point.
    pub choice: RegionChoice,
    /// `n_{m0}, n_{m1}, n_{m2}` of a vertex inserted there, by corner (from the region).
    pub counts: [i64; 3],
    /// The point on the input face of the region, as App. A computes it, before any
    /// treatment of coordinates within rounding of zero ([`RecoveredPoint`]).
    pub recovered: RecoveredPoint,
    /// The crossings of the three sides: their number and the reports of the strips
    /// (clamps, SPEC I5; mismatches).
    pub sides: StripReport,
    /// Crossing parameters along the sides (`distance / ℓ`) that fell outside `[0, 1]` by
    /// rounding and were clamped (SPEC I5), and the largest such clamp.
    pub parameter_clamps: usize,
    /// See [`parameter_clamps`](Self::parameter_clamps).
    pub max_parameter_clamp: f64,
}

/// What [`IntegerTriangulation::split_face`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceSplitReport {
    /// The new vertex, its place, the new edges and faces (as the signposts'
    /// [`Insertion`]; `snapped` is always `None`: a vertex inserted in a face stays in its
    /// input face, SPEC §2.2).
    pub insertion: Insertion,
    /// The region and the recovery of the place.
    pub location: FaceLocation,
    /// The place recovered by App. A and how it was made a point strictly inside its input
    /// face (SPEC I16).
    pub recover: RecoverReport,
}

/// One side of the face being split: its crossings in order from the side's tail.
pub(crate) struct Side {
    /// `n` of the side's edge.
    pub(crate) n: i64,
    /// The parameter of crossing `p` along the side from its tail (`distance / ℓ`).
    pub(crate) params: Vec<f64>,
    /// The place of crossing `p` on the input mesh (a point of its input edge).
    pub(crate) points: Vec<SurfacePoint>,
    /// The path of the side's edge over the input mesh (`None` along an input edge).
    pub(crate) trace: Option<EdgeTrace>,
    /// Whether the side runs against its edge's direction (its crossings are those of the
    /// edge in reverse).
    pub(crate) reversed: bool,
}

/// A point of the boundary of a region: a corner of the face, or crossing `p` of side `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RegionPoint {
    Corner(usize),
    Crossing(usize, usize),
}

/// A chord of the face: the `t`-th arc around corner `k`, or the `s`-th arc of the fan.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Chord {
    Around(usize, usize),
    Fan(usize),
}

/// The face being split or queried, with its arcs.
pub(crate) struct FaceArcs {
    pub(crate) face: usize,
    pub(crate) arrangement: FaceArrangement,
    pub(crate) layout: [[f64; 2]; 3],
    pub(crate) sides: [Side; 3],
    pub(crate) report: StripReport,
    /// Crossing parameters `distance / ℓ` outside `[0, 1]` that were clamped, and the largest
    /// clamp.
    pub(crate) parameter_clamps: usize,
    pub(crate) max_parameter_clamp: f64,
}

impl FaceArcs {
    /// Every arc of the face: those around each corner, then the fan.
    pub(crate) fn chords(&self) -> Vec<Chord> {
        let mut out: Vec<Chord> = (0..3)
            .flat_map(|k| (0..self.c(k)).map(move |t| Chord::Around(k, t)))
            .collect();
        out.extend((0..self.arrangement.fan as usize).map(Chord::Fan));
        out
    }

    /// The two ends of an arc: the `t`-th arc around corner `k` runs from crossing `t` of side
    /// `k` to crossing `n_{k+2} − 1 − t` of side `k + 2`; the `s`-th arc of the fan from the
    /// emanating corner `κ` to crossing `c_{κ+1} + s` of side `κ + 1`.
    pub(crate) fn chord_ends(&self, chord: Chord) -> (RegionPoint, RegionPoint) {
        match chord {
            Chord::Around(k, t) => {
                let k2 = (k + 2) % 3;
                (
                    RegionPoint::Crossing(k, t),
                    RegionPoint::Crossing(k2, self.n(k2) - 1 - t),
                )
            }
            Chord::Fan(s) => {
                let kappa = self.arrangement.kappa;
                let j1 = (kappa + 1) % 3;
                (
                    RegionPoint::Corner(kappa),
                    RegionPoint::Crossing(j1, self.c(j1) + s),
                )
            }
        }
    }

    /// `n⁺` of side `k` as an index bound.
    pub(crate) fn n(&self, k: usize) -> usize {
        self.sides[k].params.len()
    }

    /// `c_k` as an index.
    pub(crate) fn c(&self, k: usize) -> usize {
        self.arrangement.coords.c[k] as usize
    }

    /// The edges of the face at corner `k` in the plane: `u = P_{k+1} − P_k`, `w = P_{k+2} − P_k`.
    fn corner_frame(&self, k: usize) -> ([f64; 2], [f64; 2]) {
        let p = self.layout[k];
        let (a, b) = (self.layout[(k + 1) % 3], self.layout[(k + 2) % 3]);
        ([a[0] - p[0], a[1] - p[1]], [b[0] - p[0], b[1] - p[1]])
    }

    /// The side test of point `v` (barycentric) against `chord`, in the frame of the chord's
    /// corner: the signed distance of the point from the chord's line (positive on its
    /// positive side: the side of the corner for an arc around a corner, of corner `κ + 1` for
    /// an arc of the fan) and its distance from the chord as a segment.
    ///
    /// The sign comes from the barycentric form of the orientation, which keeps its accuracy
    /// relative to the chord's parameters however close the chord is to its corner: for the arc
    /// from the point at `a` along side `k` to the point at `b` along side `k + 2` (both from
    /// corner `k`), `f = ab − b v_{k+1} − a v_{k+2}`; for the fan arc from corner `κ` to the point
    /// at `β` along side `κ + 1`, `g = β v_{κ+1} − (1 − β) v_{κ+2}`. Both are the orientation
    /// determinant of the chord's ends and the point in barycentric coordinates, i.e. twice the
    /// signed area over twice the face's area; positions in the plane are computed relative to
    /// the corner (no cancellation with the corner's coordinates).
    fn chord_test(&self, chord: Chord, v: [f64; 3]) -> (f64, f64) {
        let (corner, rel_a, rel_b, det) = match chord {
            Chord::Around(k, t) => {
                let (k1, k2) = ((k + 1) % 3, (k + 2) % 3);
                let a = self.sides[k].params[t];
                let b = 1.0 - self.sides[k2].params[self.n(k2) - 1 - t];
                let (u, w) = self.corner_frame(k);
                (
                    k,
                    [a * u[0], a * u[1]],
                    [b * w[0], b * w[1]],
                    a * b - b * v[k1] - a * v[k2],
                )
            }
            Chord::Fan(s) => {
                let kappa = self.arrangement.kappa;
                let (k1, k2) = ((kappa + 1) % 3, (kappa + 2) % 3);
                let beta = self.sides[k1].params[self.c(k1) + s];
                let (u, w) = self.corner_frame(kappa);
                (
                    kappa,
                    [
                        (1.0 - beta) * u[0] + beta * w[0],
                        (1.0 - beta) * u[1] + beta * w[1],
                    ],
                    [0.0, 0.0],
                    beta * v[k1] - (1.0 - beta) * v[k2],
                )
            }
        };
        let (u, w) = self.corner_frame(corner);
        let twice_area = u[0] * w[1] - u[1] * w[0];
        let (c1, c2) = ((corner + 1) % 3, (corner + 2) % 3);
        let m = [v[c1] * u[0] + v[c2] * w[0], v[c1] * u[1] + v[c2] * w[1]];
        let segment = segment_distance((rel_a, rel_b), m);
        let len = (rel_b[0] - rel_a[0]).hypot(rel_b[1] - rel_a[1]);
        let signed = if len > 0.0 {
            det * twice_area / len
        } else {
            // A chord collapsed onto its corner: the point is beyond it.
            -segment
        };
        (signed, segment)
    }

    /// The chords that bound `region`, each with the sign of the region's side (`+1`: the
    /// positive side).
    pub(crate) fn bounds(&self, region: Region) -> Vec<(Chord, f64)> {
        let mut out = Vec::with_capacity(5);
        match region {
            Region::Corner { corner, index } => {
                let t = index as usize;
                out.push((Chord::Around(corner, t), 1.0));
                if t > 0 {
                    out.push((Chord::Around(corner, t - 1), -1.0));
                }
            }
            Region::Middle { .. } if self.arrangement.fan == 0 => {
                for k in 0..3 {
                    if self.c(k) > 0 {
                        out.push((Chord::Around(k, self.c(k) - 1), -1.0));
                    }
                }
            }
            Region::Middle { index } => {
                // The fan region `s` lies between the fan arcs `s − 1` and `s`; the first also
                // touches the last arc around corner `κ + 1`, the last the one around `κ + 2`
                // (corner `κ` has `c = 0`).
                let kappa = self.arrangement.kappa;
                let (s, e) = (index as usize, self.arrangement.fan as usize);
                for (j, end) in [((kappa + 1) % 3, 0), ((kappa + 2) % 3, e)] {
                    if s == end && self.c(j) > 0 {
                        out.push((Chord::Around(j, self.c(j) - 1), -1.0));
                    }
                }
                if s < e {
                    out.push((Chord::Fan(s), 1.0));
                }
                if s > 0 {
                    out.push((Chord::Fan(s - 1), -1.0));
                }
            }
        }
        out
    }

    /// The corners of the polygon `region` (3 to 6 of them, SPEC I16), in order around it.
    pub(crate) fn region_points(&self, region: Region) -> Vec<RegionPoint> {
        use RegionPoint::{Corner, Crossing};
        let mut out = Vec::with_capacity(6);
        match region {
            Region::Corner { corner: k, index } => {
                let (t, k2) = (index as usize, (k + 2) % 3);
                let n2 = self.n(k2);
                if t == 0 {
                    out.extend([Corner(k), Crossing(k, 0), Crossing(k2, n2 - 1)]);
                } else {
                    out.extend([
                        Crossing(k, t - 1),
                        Crossing(k, t),
                        Crossing(k2, n2 - 1 - t),
                        Crossing(k2, n2 - t),
                    ]);
                }
            }
            Region::Middle { index } if self.arrangement.fan == 0 => {
                let _ = index;
                for k in 0..3 {
                    let k2 = (k + 2) % 3;
                    if self.c(k) > 0 {
                        out.push(Crossing(k, self.c(k) - 1));
                        out.push(Crossing(k2, self.n(k2) - self.c(k)));
                    } else {
                        out.push(Corner(k));
                    }
                }
            }
            Region::Middle { index } => {
                let kappa = self.arrangement.kappa;
                let (j1, j2) = ((kappa + 1) % 3, (kappa + 2) % 3);
                let (s, e) = (index as usize, self.arrangement.fan as usize);
                out.push(Corner(kappa));
                if s == 0 {
                    if self.c(j1) > 0 {
                        out.push(Crossing(kappa, self.n(kappa) - self.c(j1)));
                        out.push(Crossing(j1, self.c(j1) - 1));
                    } else {
                        out.push(Corner(j1));
                    }
                }
                if s > 0 {
                    out.push(Crossing(j1, self.c(j1) + s - 1));
                }
                if s < e {
                    out.push(Crossing(j1, self.c(j1) + s));
                }
                if s == e {
                    if self.c(j2) > 0 {
                        out.push(Crossing(j1, self.n(j1) - self.c(j2)));
                        out.push(Crossing(j2, self.c(j2) - 1));
                    } else {
                        out.push(Corner(j2));
                    }
                }
            }
        }
        out
    }

    /// A gap of `region` on the face's boundary: `(side, slot)`, the piece of side `side`
    /// between its crossings `slot − 1` and `slot` (from the side's tail), on a side that does
    /// not run along an input edge; `None` when every side runs along one (the face is an
    /// input face).
    pub(crate) fn region_gap(&self, region: Region) -> Option<(usize, usize)> {
        match region {
            Region::Corner { corner, index } => Some((corner, index as usize)),
            Region::Middle { index } if self.arrangement.fan > 0 => {
                let j1 = (self.arrangement.kappa + 1) % 3;
                Some((j1, self.c(j1) + index as usize))
            }
            Region::Middle { .. } => (0..3)
                .find(|&k| self.sides[k].n >= 0)
                .map(|k| (k, self.c(k))),
        }
    }

    /// The barycentric coordinates on the face of a region point.
    pub(crate) fn bary_of(&self, p: RegionPoint) -> [f64; 3] {
        let mut b = [0.0; 3];
        match p {
            RegionPoint::Corner(k) => b[k] = 1.0,
            RegionPoint::Crossing(k, i) => {
                let s = self.sides[k].params[i];
                b[k] = 1.0 - s;
                b[(k + 1) % 3] = s;
            }
        }
        b
    }

    /// The region of the point `v` (barycentric), by the largest margin (I15), and the distance
    /// to the nearest arc (I6).
    pub(crate) fn choose_region(&self, v: [f64; 3]) -> RegionChoice {
        let a = &self.arrangement;
        let l = [0, 1, 2].map(|k| {
            let (p, q) = (self.layout[k], self.layout[(k + 1) % 3]);
            (q[0] - p[0]).hypot(q[1] - p[1])
        });
        let unit = l[0].min(l[1]).min(l[2]);
        // Side tests against every chord: (signed distance, segment distance).
        let around: [Vec<(f64, f64)>; 3] = [0, 1, 2].map(|k| {
            (0..self.c(k))
                .map(|t| self.chord_test(Chord::Around(k, t), v))
                .collect()
        });
        let fan: Vec<(f64, f64)> = (0..a.fan as usize)
            .map(|s| self.chord_test(Chord::Fan(s), v))
            .collect();
        let test = |chord: Chord| match chord {
            Chord::Around(k, t) => around[k][t],
            Chord::Fan(s) => fan[s],
        };
        let margin = |region: Region| {
            self.bounds(region)
                .into_iter()
                .map(|(chord, sign)| sign * test(chord).0)
                .fold(
                    f64::INFINITY,
                    |x, y| if y < x || y.is_nan() { y } else { x },
                )
        };
        let key = |x: f64| if x.is_nan() { f64::NEG_INFINITY } else { x };
        let mut best: Option<(Region, f64)> = None;
        let mut second: Option<(Region, f64)> = None;
        for region in a.regions() {
            let g = key(margin(region));
            match best {
                Some((_, b)) if g <= b => {
                    if second.is_none_or(|(_, s)| g > s) {
                        second = Some((region, g));
                    }
                }
                _ => {
                    second = best;
                    best = Some((region, g));
                }
            }
        }
        let (region, g) = best.expect("every face has a middle region");
        let bounds = self.bounds(region);
        let curve_distance = (!bounds.is_empty()).then(|| {
            bounds
                .iter()
                .map(|&(chord, _)| test(chord).1)
                .fold(f64::INFINITY, f64::min)
                / unit
        });
        RegionChoice {
            region,
            margin: g / unit,
            runner_up: second.map(|(r, s)| (r, s / unit)),
            curve_distance,
            n_regions: a.n_regions(),
        }
    }
}

/// The distance of `m` from the segment `(a, b)`.
fn segment_distance((a, b): ([f64; 2], [f64; 2]), m: [f64; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let w = [m[0] - a[0], m[1] - a[1]];
    let dd = d[0] * d[0] + d[1] * d[1];
    let s = if dd > 0.0 {
        ((w[0] * d[0] + w[1] * d[1]) / dd).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (w[0] - s * d[0]).hypot(w[1] - s * d[1])
}

impl IntegerTriangulation {
    /// The arcs of face `f`: its arrangement from the integers, its layout in the plane, and
    /// the crossings of its sides from `ExtractIntrinsicEdge`.
    pub(crate) fn face_arcs(&self, f: usize) -> Result<FaceArcs, IntegerError> {
        let nf = self.complex.n_faces();
        if f >= nf {
            return Err(IntegerError::IndexOutOfRange {
                what: "face",
                index: f,
                bound: nf,
            });
        }
        let n = self.face_normal(f);
        let arrangement = FaceArrangement::of(n).ok_or_else(|| {
            IntegerError::invalid(
                ValidityItem::FaceParity,
                Place::Face(f as u32),
                format!("face {f} with side counts {n:?} is not valid"),
            )
        })?;
        let l = self.face_lengths(f);
        let layout =
            triangle_layout(l).ok_or(IntegerError::Intrinsic(IntrinsicError::DegenerateFace {
                surface: "intrinsic triangulation",
                face: f as u32,
                lengths: l,
            }))?;
        let mut report = StripReport::default();
        let mut clamps = (0usize, 0.0f64);
        let mut side = |k: usize| -> Result<Side, IntegerError> {
            let h = 3 * f + k;
            let e = self.complex.edge(h);
            let reversed = self.complex.edge_halfedge(e) != h;
            let n = self.normal[e];
            if n < 0 {
                return Ok(Side {
                    n,
                    params: Vec::new(),
                    points: Vec::new(),
                    trace: None,
                    reversed,
                });
            }
            let (trace, r) =
                self.extract_intrinsic_edge_with_report(e, &EdgeTraceOptions::default())?;
            report.merge(&r);
            let x = &trace.crossings;
            if x.len() != n as usize {
                return Err(IntegerError::invalid(
                    ValidityItem::Extraction,
                    Place::Edge(e as u32),
                    format!("{} crossings extracted for n = {n}", x.len()),
                ));
            }
            let len = self.lengths[e];
            let order = |p: usize| if reversed { x.len() - 1 - p } else { p };
            // `distance / ℓ`, clamped into [0, 1] (the strip's length differs from ℓ by
            // rounding; SPEC I5) and counted.
            let params = (0..x.len())
                .map(|p| {
                    let s = x[order(p)].distance / len;
                    let s = if reversed { 1.0 - s } else { s };
                    let c = s.clamp(0.0, 1.0);
                    if c != s {
                        clamps.0 += 1;
                        clamps.1 = clamps.1.max((c - s).abs());
                    }
                    c
                })
                .collect();
            let points = (0..x.len())
                .map(|p| SurfacePoint::Edge {
                    edge: x[order(p)].edge,
                    t: x[order(p)].t,
                })
                .collect();
            Ok(Side {
                n,
                params,
                points,
                trace: Some(trace),
                reversed,
            })
        };
        let sides = [side(0)?, side(1)?, side(2)?];
        Ok(FaceArcs {
            face: f,
            arrangement,
            layout,
            sides,
            report,
            parameter_clamps: clamps.0,
            max_parameter_clamp: clamps.1,
        })
    }

    /// The point query of SPEC §4.9 inside face `f` of the triangulation, without inserting
    /// anything: the arcs of the face, the region of the point `bary` (I15), the counts a
    /// vertex inserted there would get, and its place on the input mesh by App. A
    /// ([`crate::recover`], every corner of the region, I16), as computed (the coordinates may
    /// be within rounding of zero, see [`RecoveredPoint`]).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] with
    /// [`IntrinsicError::InvalidPoint`] for coordinates that are not positive and summing to
    /// one, or a degenerate face; the errors of
    /// [`extract_intrinsic_edge`](Self::extract_intrinsic_edge); [`IntegerError::RecoverFailed`]
    /// when App. A has no solution or a corner of the region is not on the input face of the
    /// region.
    pub fn locate_in_face(&self, f: usize, bary: [f64; 3]) -> Result<FaceLocation, IntegerError> {
        let nf = self.complex.n_faces();
        if f >= nf {
            return Err(IntegerError::IndexOutOfRange {
                what: "face",
                index: f,
                bound: nf,
            });
        }
        let at = SurfacePoint::Face {
            face: f as u32,
            bary,
        };
        if !bary.iter().all(|&b| b.is_finite() && b > 0.0)
            || ((bary[0] + bary[1] + bary[2]) - 1.0).abs() > 1e-9
        {
            return Err(IntegerError::Intrinsic(IntrinsicError::InvalidPoint {
                point: at,
                reason: "the barycentric coordinates of a point inside a face must be positive \
                         and sum to one",
            }));
        }
        let arcs = self.face_arcs(f)?;
        let choice = arcs.choose_region(bary);
        let counts =
            arcs.arrangement
                .split_counts(choice.region)
                .ok_or(IntegerError::Overflow {
                    edge: self.complex.edge(3 * f) as u32,
                })?;
        let recovered = recover_in_input_face(self, &arcs, choice.region, bary)?;
        Ok(FaceLocation {
            face: f as u32,
            arrangement: arcs.arrangement,
            choice,
            counts,
            recovered,
            sides: arcs.report,
            parameter_clamps: arcs.parameter_clamps,
            max_parameter_clamp: arcs.max_parameter_clamp,
        })
    }

    /// Inserts a vertex at the point `bary` of face `f` (paper Sec. 3.3, Supp. Alg. 5, 6;
    /// SPEC §4.4) and joins it to the three corners, with the default options (I6 never
    /// refuses). Ids as [`Complex::split_face`](ddg_intrinsic::Complex::split_face): the face
    /// keeps its id and becomes `(a, b, m)`, two faces, three edges (from corners `a`, `b`, `c`
    /// to `m`) and one vertex are appended. See [`split_face`](Self::split_face).
    ///
    /// # Errors
    /// As [`split_face`](Self::split_face).
    pub fn insert_vertex(&mut self, f: usize, bary: [f64; 3]) -> Result<Insertion, IntegerError> {
        self.split_face(f, bary, &FaceSplitOptions::default())
            .map(|r| r.insertion)
    }

    /// `SplitFace` (Supp. Alg. 5, 6; SPEC §4.4, I6, I15, I16, I25): see the
    /// [module documentation](crate::insert). The new edges get the lengths of
    /// [`plan_face_split`] (bitwise as the signposts), the normal coordinates of the region
    /// the point is in, the roundabouts of Eq. 4's analogue at input-vertex corners; the new
    /// vertex gets the place of App. A, a point strictly inside its input face.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] with the refusals of [`plan_face_split`] (index, point,
    /// degenerate face, [`InsertRefusal::TooThin`](ddg_intrinsic::InsertRefusal::TooThin)) and
    /// [`IntrinsicError::InvalidOptions`] for a `too_close` that is not finite and not negative;
    /// [`IntegerError::TooCloseToCurve`] (I6); [`IntegerError::RecoverFailed`] when App. A
    /// gives a point outside its input face by more than rounding (I16); the errors of
    /// [`locate_in_face`](Self::locate_in_face); [`IntegerError::InvalidCoordinates`] for a
    /// roundabout that gives none for a new edge. Nothing is changed in every case.
    pub fn split_face(
        &mut self,
        f: usize,
        bary: [f64; 3],
        options: &FaceSplitOptions,
    ) -> Result<FaceSplitReport, IntegerError> {
        if !(options.too_close.is_finite() && options.too_close >= 0.0) {
            return Err(IntegerError::Intrinsic(IntrinsicError::InvalidOptions {
                reason: format!(
                    "FaceSplitOptions::too_close must be finite and not negative, got {}",
                    options.too_close
                ),
            }));
        }
        let plan = plan_face_split(&self.complex, &self.lengths, f, bary)?;
        let location = self.locate_in_face(f, bary)?;
        if let Some(d) = location.choice.curve_distance
            && d < options.too_close
        {
            return Err(IntegerError::TooCloseToCurve {
                face: f as u32,
                distance: d,
                limit: options.too_close,
            });
        }
        let (place, recover) = location.recovered.inside_input_face(f)?;
        // The roundabouts of the new halfedges from the corners, before anything changes.
        let c = &self.complex;
        let mut roundabouts = [NO_ROUNDABOUT; 3];
        for (k, slot) in roundabouts.iter_mut().enumerate() {
            let h = 3 * f + k;
            let a = c.from_vertex(h);
            if !self.is_input_vertex(a) {
                continue;
            }
            let deg = self.local.degree(a);
            let r = self.roundabouts[h];
            let before = location
                .arrangement
                .emanating_before(k, location.choice.region);
            let x = (r < deg)
                .then(|| {
                    roundabout_update(
                        r as i64,
                        before,
                        n_minus(self.normal[c.edge(h)]),
                        self.local.modulus(a),
                    )
                })
                .flatten();
            *slot = match x {
                Some(x) if (0..deg as i64).contains(&x) => x as u32,
                _ => {
                    return Err(IntegerError::invalid(
                        ValidityItem::Ranges,
                        Place::Halfedge(h as u32),
                        format!(
                            "the roundabout {r} of halfedge {h} at input vertex {a} (deg₀ {deg}) \
                             gives no roundabout for the new edge to the inserted vertex"
                        ),
                    ));
                }
            };
        }
        // Change the state.
        let split = self
            .complex
            .split_face(f)
            .expect("plan_face_split checked the face");
        split.move_halfedge_data(&mut self.roundabouts, NO_ROUNDABOUT);
        for (k, &e) in split.new_edges.iter().enumerate() {
            self.lengths.push(plan.distances[k]);
            self.normal.push(location.counts[k]);
            debug_assert_eq!(self.lengths.len(), e as usize + 1);
            let to_m = self.complex.edge_halfedge(e as usize);
            self.roundabouts[to_m] = roundabouts[k];
        }
        for h in split.outgoing {
            self.roundabouts[h as usize] = NO_ROUNDABOUT;
        }
        self.locations.push(place);
        let insertion = Insertion {
            vertex: split.vertex,
            location: place,
            new_edges: split.new_edges.to_vec(),
            new_faces: split.new_faces.to_vec(),
            halves: None,
            snapped: None,
        };
        Ok(FaceSplitReport {
            insertion,
            location,
            recover,
        })
    }
}
