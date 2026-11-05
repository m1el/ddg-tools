/-!
# Integer coordinates: the integer formulas (`crates/ddg-integer-coords/SPEC.md`)

"Integer Coordinates for Intrinsic Geometry Processing", Gillespie, Sharp, Crane, ACM TOG 40(6), 2021
("the paper"; Supp. = its supplement). Written from the paper's text (Sec. 2.3–3.7, Supp. B), not from
the authors' code. Everything here is exact integer arithmetic (`Int`, `Nat`); there is no rounding.
`DDG/NormalCoords/*.lean` proves theorems about these very definitions.

**Face convention** (the halfedge convention of the workspace, `h = 3f + k` runs from corner `k` to
corner `k + 1`). A face has corners `0, 1, 2` and sides `0, 1, 2`; side `k` runs from corner `k` to
corner `k + 1`, and its crossings are numbered `0, 1, …` from corner `k` (the paper's `p`-th
crossing along the halfedge from its tail). Corner `k` lies between side `k + 2` (which ends at it)
and side `k` (which starts at it); the side opposite corner `k` is side `k + 1`. Indices are taken
modulo 3 (`next`, `prev`).

**Never `/` on a negative integer** (SPEC I14): the numerator of Eq. 1 is never negative
(`DDG.IntegerCoords.cornerNumerator_nonneg`), and `cornerCoord?` refuses an odd one instead
of rounding it.
-/

namespace Ddg.IntegerCoords

/-! ## Signs (Sec. 2.3) -/

/-- `n⁺ = max(n, 0)`: the number of transversal crossings. [formula: ic-normal-sign] -/
def nPlus (n : Int) : Int := max n 0

/-- `n⁻ = −min(n, 0)`: `1` on an edge that runs along an input edge (`n = −1`), else `0`.
[formula: ic-normal-sign] -/
def nMinus (n : Int) : Int := -(min n 0)

/-! ## Corner coordinates and emanating counts (Eq. 1, 2) -/

/-- Eq. 2, `e^{ij}_k = max(0, n⁺_ij − n⁺_jk − n⁺_ki)`: arcs from corner `k` to the opposite side
`ij`. Arguments: the opposite side, then the two sides at the corner (in either order).
[formula: ic-emanating] -/
def emanating (opp a b : Int) : Int := max 0 (nPlus opp - nPlus a - nPlus b)

/-- Twice Eq. 1 at corner `k` of `ijk`: `max(0, n⁺_jk + n⁺_ki − n⁺_ij) − e^{jk}_i − e^{ki}_j`.
Arguments: the opposite side `ij`, then the two sides `jk`, `ki` at the corner.
[formula: ic-corner-coordinates] -/
def cornerNumerator (opp a b : Int) : Int :=
  max 0 (nPlus a + nPlus b - nPlus opp) - emanating a opp b - emanating b opp a

/-- Eq. 1, `c^{ij}_k`: the numerator halved. Only meaningful when the numerator is even
(`faceValid`); the numerator is never negative, so this `/` never sees a negative integer.
[formula: ic-corner-coordinates] -/
def cornerCoord (opp a b : Int) : Int := cornerNumerator opp a b / 2

/-- Eq. 1 checked: `none` when the numerator is negative (impossible) or odd (an invalid face).
[formula: ic-corner-coordinates] -/
def cornerCoord? (opp a b : Int) : Option Int :=
  let m := cornerNumerator opp a b
  if m < 0 then none else if m % 2 != 0 then none else some (m / 2)

/-- A face `(n₀, n₁, n₂)` (side `k` from corner `k` to `k + 1`) is valid (SPEC §2.4 item 1) when
the three numerators of Eq. 1 are even. [formula: ic-face-validity] -/
def faceValid (n0 n1 n2 : Int) : Bool :=
  cornerNumerator n1 n0 n2 % 2 == 0 && cornerNumerator n2 n1 n0 % 2 == 0 &&
    cornerNumerator n0 n2 n1 % 2 == 0

/-- `k + 1 (mod 3)`. -/
def next (k : Nat) : Nat := (k + 1) % 3
/-- `k + 2 = k − 1 (mod 3)`. -/
def prev (k : Nat) : Nat := (k + 2) % 3

/-- Side counts of a face as a function of the side index (`0, 1, 2`; anything else is side 2). -/
def side (n0 n1 n2 : Int) (k : Nat) : Int :=
  if k % 3 == 0 then n0 else if k % 3 == 1 then n1 else n2

/-- Corner coordinate of corner `k` of the face `(n₀, n₁, n₂)`: opposite side `k + 1`, sides `k`
and `k + 2` at the corner. [formula: ic-corner-coordinates] -/
def faceCorner (n0 n1 n2 : Int) (k : Nat) : Int :=
  cornerCoord (side n0 n1 n2 (next k)) (side n0 n1 n2 k) (side n0 n1 n2 (prev k))

/-- Emanating count of corner `k` (arcs from corner `k` to side `k + 1`). [formula: ic-emanating] -/
def faceEmanating (n0 n1 n2 : Int) (k : Nat) : Int :=
  emanating (side n0 n1 n2 (next k)) (side n0 n1 n2 k) (side n0 n1 n2 (prev k))

/-! ## Crossings (Sec. 2.5) -/

/-- The reversal of a combinatorial crossing: `(ij→, p) ↦ (ji→, n_ij − p − 1)`; the position part.
[formula: ic-reverse-crossing] -/
def reverseCrossing (n p : Int) : Int := n - p - 1

/-! ## Edge flip (Eq. 3, 4) -/

/-- Absolute value of an integer. -/
def iabs (x : Int) : Int := if x < 0 then -x else x

/-- Twice Eq. 3 without `n⁻_ij`. Edge `ij` with the faces `ijk` and `jil`; arguments
`n_ij, n_jk, n_ki, n_il, n_lj`. Every corner quantity is that of one of the two old faces:
`c^{ji}_l` is corner `l` of `jil` (the paper prints `c^{jk}_l`, a misprint: `l` and `jk` share no
face, and the reference's code uses corner `l` of `jil`), `c^{ij}_k` corner `k` of `ijk`,
`c^{il}_j`, `e^{il}_j` corner `j` of `jil` and `c^{ki}_j`, `e^{ki}_j` corner `j` of `ijk`, and so on.
[formula: ic-flip-normal] -/
def flipDoubled (nij njk nki nil nlj : Int) : Int :=
  let cl := cornerCoord nij nil nlj          -- c^{ji}_l, corner l of jil
  let ck := cornerCoord nij njk nki          -- c^{ij}_k, corner k of ijk
  let cjL := cornerCoord nil nij nlj         -- c^{il}_j, corner j of jil
  let cjR := cornerCoord nki nij njk         -- c^{ki}_j, corner j of ijk
  let ciL := cornerCoord nlj nij nil         -- c^{lj}_i, corner i of jil
  let ciR := cornerCoord njk nij nki         -- c^{jk}_i, corner i of ijk
  let el := emanating nij nil nlj            -- e^{ji}_l
  let ek := emanating nij njk nki            -- e^{ij}_k
  let eiL := emanating nlj nij nil           -- e^{lj}_i
  let eiR := emanating njk nij nki           -- e^{jk}_i
  let ejL := emanating nil nij nlj           -- e^{il}_j
  let ejR := emanating nki nij njk           -- e^{ki}_j
  2 * cl + 2 * ck + iabs (cjL - cjR) + iabs (ciL - ciR) - el - ek +
    2 * (eiL + eiR + ejL + ejR)

/-- Eq. 3: the normal coordinate of the new edge `kl` after flipping `ij` (faces `ijk`, `jil`).
[formula: ic-flip-normal] -/
def flipNormal (nij njk nki nil nlj : Int) : Int :=
  flipDoubled nij njk nki nil nlj / 2 + nMinus nij

/-- Eq. 4: `r_kl = (r_ki + e^{il}_k + n⁻_ki) mod deg₀(k)`; `e^{il}_k` is the emanating count of
corner `k` in the new face `kil`. With `deg = 0` (a boundary vertex, whose numbering does not wrap,
SPEC I3) no reduction. [formula: ic-roundabout-update] -/
def roundaboutUpdate (r eNew nMinusKi : Int) (deg : Nat) : Int :=
  if deg == 0 then r + eNew + nMinusKi else (r + eNew + nMinusKi) % deg

/-! ## Edge split along a shared edge (Eq. 6) -/

/-- Eq. 6: `n_mk = max(n_ki, n_jk, 0)` for the new edge from the split point `m` of an edge `ij`
with `n_ij = −1` to the opposite corner `k` (the halves keep `n = −1`). [formula: ic-split-edge-normal] -/
def splitEdgeNormal (nki njk : Int) : Int := max (max nki njk) 0

/-! ## Tracing one step (Supp. Alg. 1) -/

/-- Where a curve that crosses an edge goes next in the face it enters. -/
inductive Step where
  /-- leaves through `ik→` at this crossing index (Fig. 7 case 3, "goes right") -/
  | toIK (p : Int)
  /-- leaves through `kj→` at this crossing index (Fig. 7 case 1, "goes left") -/
  | toKJ (p : Int)
  /-- ends at the opposite vertex `k` (Fig. 7 case 2) -/
  | vertex
  deriving DecidableEq, Repr, Inhabited

/-- One step of Alg. 1: the curve crosses `ij→` at index `p` (from `i`) and enters the face `jik`
of the twin, with opposite vertex `k`. `c^{kj}_i` (corner `i` of `jik`) and `c^{ik}_j` (corner `j`).
[formula: ic-trace-step] -/
def traceStep (nij nik nkj p : Int) : Step :=
  if p < cornerCoord nkj nij nik then .toIK p
  else if nPlus nij - cornerCoord nik nij nkj ≤ p then .toKJ (nPlus nkj + p - nPlus nij)
  else .vertex

/-! ## Identifying a traced edge (Supp. Alg. 3, the integer part) -/

/-- Alg. 3 line 9: the curve leaves `a ∈ V⁰` and first crosses `ij→` at `p`, with `ija` the face of
`ij→` and `ai→ = Next(Next(ij→))` (the paper calls it `aj`). Its local index around `a` is
`r_ai + n⁻_ai + (p − c^{ja}_i)`, reduced modulo `deg₀(a)` at an interior vertex (`deg = 0`: a
boundary vertex, no reduction). Arguments `r_ai, n_ai, n_ij, n_ja, p, deg`.
[formula: ic-identify-edge] -/
def identifyIndex (rAi nAi nIj nJa p : Int) (deg : Nat) : Int :=
  let k := rAi + nMinus nAi + (p - cornerCoord nJa nIj nAi)
  if deg == 0 then k else k % deg

/-! ## The arcs of one face (Supp. Alg. 9, Fig. 11) -/

/-- An endpoint of an arc inside a face. -/
inductive Pt where
  | corner (k : Nat)
  | side (k : Nat) (p : Nat)
  deriving DecidableEq, Repr, Inhabited

/-- The other end of the arc through crossing `p` of side `k` (`0 ≤ p < n_k`), for the face with
side counts `n` and corner coordinates `c` (`n`, `c` as the positive parts): the first `c_k` crossings
(from corner `k`) cut corner `k` and reach side `k + 2` in reverse order, the last `c_{k+1}` cut
corner `k + 1` and reach side `k + 1` in reverse order, and the ones between come from the opposite
corner `k + 2`. This is the pairing of Alg. 9 in both of its cases (case 2 is the one where the
middle block is not empty); it is the only valid one (`DDG.IntegerCoords.FaceModel.eq_canon`).
[formula: ic-face-decoding] -/
def partner (n c : Nat → Nat) (k p : Nat) : Pt :=
  if p < c k then .side (prev k) (n (prev k) - 1 - p)
  else if p + c (next k) < n k then .corner (prev k)
  else .side (next k) (n k - 1 - p)

/-! ## Face split: regions and counts (Supp. Alg. 5, 6; SPEC I15) -/

/-- A region of the arrangement of a face's arcs. `corner k t` (`t < c_k`) is the region between the
`t − 1`-st and the `t`-th arc around corner `k` (counted from the corner; `t = 0` contains the
corner). `middle s` is what is left: with no emanating arcs the one central region (`s = 0`); with
`e > 0` arcs from corner `κ`, the `s`-th region of the fan (`0 ≤ s ≤ e`), counted from corner
`κ + 1`. -/
inductive Region where
  | corner (k t : Nat)
  | middle (s : Nat)
  deriving DecidableEq, Repr, Inhabited

/-- The decoded face: `n⁺` per side, `c` per corner, and the emanating corner with its count
(`κ = 0` and `e = 0` when nothing emanates). `none` for an invalid face. -/
structure Decoded where
  n : Nat → Nat
  c : Nat → Nat
  kappa : Nat
  e : Nat

/-- Decode a face `(n₀, n₁, n₂)` (Eq. 1, 2). [formula: ic-face-decoding] -/
def decode? (n0 n1 n2 : Int) : Option Decoded :=
  if !faceValid n0 n1 n2 then none else
  let nn : Nat → Nat := fun k => (nPlus (side n0 n1 n2 k)).toNat
  let cc : Nat → Nat := fun k => (faceCorner n0 n1 n2 k).toNat
  let es := [0, 1, 2].map fun k => (k, (faceEmanating n0 n1 n2 k).toNat)
  match es.find? (fun ke => ke.2 > 0) with
  | some (k, e) => some { n := nn, c := cc, kappa := k, e := e }
  | none => some { n := nn, c := cc, kappa := 0, e := 0 }

/-- All regions of a decoded face, in a fixed order: the corner regions of corners `0, 1, 2`, then
the middle regions. [formula: ic-face-split-normal] -/
def regions (d : Decoded) : List Region :=
  ([0, 1, 2].flatMap fun k => (List.range (d.c k)).map fun t => Region.corner k t) ++
    (List.range (d.e + 1)).map Region.middle

/-- Is a region one of the face's? -/
def Region.valid (d : Decoded) : Region → Bool
  | .corner k t => k < 3 && t < d.c k
  | .middle s => s ≤ d.e

/-- The normal coordinates `n_{m0}, n_{m1}, n_{m2}` of the edges from a vertex `m` inserted in
region `R` to the corners: the number of arcs that separate `R` from each corner (an arc that ends at
the corner does not). Corner region `(k, t)`: `t` arcs to corner `k`, and `c_j + c_k − t` to the
other corners `j`, plus the `e` fan arcs when `{k, j} = {κ + 1, κ + 2}`. Middle region `s`:
`c_j`, plus `s` (for `j = κ + 1`) or `e − s` (for `j = κ + 2`) fan arcs.
[formula: ic-face-split-normal] -/
def splitCount (d : Decoded) (R : Region) (j : Nat) : Nat :=
  match R with
  | .corner k t =>
    if j == k then t
    else
      let fan := if (k == next d.kappa && j == prev d.kappa) ||
          (k == prev d.kappa && j == next d.kappa) then d.e else 0
      d.c j + d.c k - t + fan
  | .middle s =>
    if j == next d.kappa then d.c j + s
    else if j == prev d.kappa then d.c j + (d.e - s)
    else d.c j

/-- Face split (I15): the three new normal coordinates for the region `R` of the face
`(n₀, n₁, n₂)`, or `none` for an invalid face or a region it does not have.
[formula: ic-face-split-normal] -/
def faceSplit? (n0 n1 n2 : Int) (R : Region) : Option (Nat × Nat × Nat) :=
  match decode? n0 n1 n2 with
  | none => none
  | some d =>
    if R.valid d then some (splitCount d R 0, splitCount d R 1, splitCount d R 2) else none

end Ddg.IntegerCoords
