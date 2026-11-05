import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.Group.Int

/-!
# Integer coordinates: the integer formulas (`crates/ddg-integer-coords/SPEC.md`)

"Integer Coordinates for Intrinsic Geometry Processing", Gillespie, Sharp, Crane, ACM TOG 40(6), 2021
("the paper"; Supp. = its supplement). Written from the paper's text (Sec. 2.3–3.7, Supp. B), not from
the authors' code. Everything here is exact integer arithmetic (`Int`, `Nat`); there is no rounding.
This is the library's native copy (namespace `DDG.NormalCoords.Formulas`) of the oracle's
`DdgOracle/NormalCoords/Formulas.lean`; the app bridge `DdgProofs/Oracle/NormalCoords.lean` proves them equal by `rfl`.

**Face convention** (the halfedge convention of the workspace, `h = 3f + k` runs from corner `k` to
corner `k + 1`). A face has corners `0, 1, 2` and sides `0, 1, 2`; side `k` runs from corner `k` to
corner `k + 1`, and its crossings are numbered `0, 1, …` from corner `k` (the paper's `p`-th
crossing along the halfedge from its tail). Corner `k` lies between side `k + 2` (which ends at it)
and side `k` (which starts at it); the side opposite corner `k` is side `k + 1`. Indices are taken
modulo 3 (`next`, `prev`).

**Never `/` on a negative integer** (SPEC I14): the numerator of Eq. 1 is never negative
(`DDG.NormalCoords.cornerNumerator_nonneg`), and `cornerCoord?` refuses an odd one instead
of rounding it.
-/

namespace DDG.NormalCoords.Formulas

/-! ## Signs (Sec. 2.3) -/

/-- `n⁺ = max(n, 0)`: the number of transversal crossings. -/
def nPlus (n : Int) : Int := max n 0

/-- `n⁻ = −min(n, 0)`: `1` on an edge that runs along an input edge (`n = −1`), else `0`. -/
def nMinus (n : Int) : Int := -(min n 0)

/-! ## Corner coordinates and emanating counts (Eq. 1, 2) -/

/-- Eq. 2, `e^{ij}_k = max(0, n⁺_ij − n⁺_jk − n⁺_ki)`: arcs from corner `k` to the opposite side
`ij`. Arguments: the opposite side, then the two sides at the corner (in either order). -/
def emanating (opp a b : Int) : Int := max 0 (nPlus opp - nPlus a - nPlus b)

/-- Twice Eq. 1 at corner `k` of `ijk`: `max(0, n⁺_jk + n⁺_ki − n⁺_ij) − e^{jk}_i − e^{ki}_j`.
Arguments: the opposite side `ij`, then the two sides `jk`, `ki` at the corner. -/
def cornerNumerator (opp a b : Int) : Int :=
  max 0 (nPlus a + nPlus b - nPlus opp) - emanating a opp b - emanating b opp a

/-- Eq. 1, `c^{ij}_k`: the numerator halved. Only meaningful when the numerator is even
(`faceValid`); the numerator is never negative, so this `/` never sees a negative integer. -/
def cornerCoord (opp a b : Int) : Int := cornerNumerator opp a b / 2

/-- Eq. 1 checked: `none` when the numerator is negative (impossible) or odd (an invalid face). -/
def cornerCoord? (opp a b : Int) : Option Int :=
  let m := cornerNumerator opp a b
  if m < 0 then none else if m % 2 != 0 then none else some (m / 2)

/-- A face `(n₀, n₁, n₂)` (side `k` from corner `k` to `k + 1`) is valid (SPEC §2.4 item 1) when
the three numerators of Eq. 1 are even. -/
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
and `k + 2` at the corner. -/
def faceCorner (n0 n1 n2 : Int) (k : Nat) : Int :=
  cornerCoord (side n0 n1 n2 (next k)) (side n0 n1 n2 k) (side n0 n1 n2 (prev k))

/-- Emanating count of corner `k` (arcs from corner `k` to side `k + 1`). -/
def faceEmanating (n0 n1 n2 : Int) (k : Nat) : Int :=
  emanating (side n0 n1 n2 (next k)) (side n0 n1 n2 k) (side n0 n1 n2 (prev k))

/-! ## Crossings (Sec. 2.5) -/

/-- The reversal of a combinatorial crossing: `(ij→, p) ↦ (ji→, n_ij − p − 1)`; the position part. -/
def reverseCrossing (n p : Int) : Int := n - p - 1

/-! ## Edge flip (Eq. 3, 4) -/

/-- Absolute value of an integer. -/
def iabs (x : Int) : Int := if x < 0 then -x else x

/-- Twice Eq. 3 without `n⁻_ij`. Edge `ij` with the faces `ijk` and `jil`; arguments
`n_ij, n_jk, n_ki, n_il, n_lj`. Every corner quantity is that of one of the two old faces:
`c^{ji}_l` is corner `l` of `jil` (the paper prints `c^{jk}_l`, a misprint: `l` and `jk` share no
face, and the reference's code uses corner `l` of `jil`), `c^{ij}_k` corner `k` of `ijk`,
`c^{il}_j`, `e^{il}_j` corner `j` of `jil` and `c^{ki}_j`, `e^{ki}_j` corner `j` of `ijk`, and so on. -/
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

/-- Eq. 3: the normal coordinate of the new edge `kl` after flipping `ij` (faces `ijk`, `jil`). -/
def flipNormal (nij njk nki nil nlj : Int) : Int :=
  flipDoubled nij njk nki nil nlj / 2 + nMinus nij

/-- Eq. 4: `r_kl = (r_ki + e^{il}_k + n⁻_ki) mod deg₀(k)`; `e^{il}_k` is the emanating count of
corner `k` in the new face `kil`. With `deg = 0` (a boundary vertex, whose numbering does not wrap,
SPEC I3) no reduction. -/
def roundaboutUpdate (r eNew nMinusKi : Int) (deg : Nat) : Int :=
  if deg == 0 then r + eNew + nMinusKi else (r + eNew + nMinusKi) % deg

/-! ## Edge split along a shared edge (Eq. 6) -/

/-- Eq. 6: `n_mk = max(n_ki, n_jk, 0)` for the new edge from the split point `m` of an edge `ij`
with `n_ij = −1` to the opposite corner `k` (the halves keep `n = −1`). -/
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
of the twin, with opposite vertex `k`. `c^{kj}_i` (corner `i` of `jik`) and `c^{ik}_j` (corner `j`). -/
def traceStep (nij nik nkj p : Int) : Step :=
  if p < cornerCoord nkj nij nik then .toIK p
  else if nPlus nij - cornerCoord nik nij nkj ≤ p then .toKJ (nPlus nkj + p - nPlus nij)
  else .vertex

/-! ## Identifying a traced edge (Supp. Alg. 3, the integer part) -/

/-- Alg. 3 line 9: the curve leaves `a ∈ V⁰` and first crosses `ij→` at `p`, with `ija` the face of
`ij→` and `ai→ = Next(Next(ij→))` (the paper calls it `aj`). Its local index around `a` is
`r_ai + n⁻_ai + (p − c^{ja}_i)`, reduced modulo `deg₀(a)` at an interior vertex (`deg = 0`: a
boundary vertex, no reduction). Arguments `r_ai, n_ai, n_ij, n_ja, p, deg`. -/
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
middle block is not empty); it is the only valid one (`DDG.NormalCoords.FaceModel.eq_canon`). -/
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

/-- Decode a face `(n₀, n₁, n₂)` (Eq. 1, 2). -/
def decode? (n0 n1 n2 : Int) : Option Decoded :=
  if !faceValid n0 n1 n2 then none else
  let nn : Nat → Nat := fun k => (nPlus (side n0 n1 n2 k)).toNat
  let cc : Nat → Nat := fun k => (faceCorner n0 n1 n2 k).toNat
  let es := [0, 1, 2].map fun k => (k, (faceEmanating n0 n1 n2 k).toNat)
  match es.find? (fun ke => ke.2 > 0) with
  | some (k, e) => some { n := nn, c := cc, kappa := k, e := e }
  | none => some { n := nn, c := cc, kappa := 0, e := 0 }

/-- All regions of a decoded face, in a fixed order: the corner regions of corners `0, 1, 2`, then
the middle regions. -/
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
`c_j`, plus `s` (for `j = κ + 1`) or `e − s` (for `j = κ + 2`) fan arcs. -/
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
`(n₀, n₁, n₂)`, or `none` for an invalid face or a region it does not have. -/
def faceSplit? (n0 n1 n2 : Int) (R : Region) : Option (Nat × Nat × Nat) :=
  match decode? n0 n1 n2 with
  | none => none
  | some d =>
    if R.valid d then some (splitCount d R 0, splitCount d R 1, splitCount d R 2) else none

end DDG.NormalCoords.Formulas

/-!
# Integer coordinates: arithmetic of the formulas (tier A)

`crates/ddg-integer-coords/SPEC.md` §2.2–2.4, §4.2, §4.5, §9.2; the paper Sec. 2.3, 2.5, 3.2, 3.4.
All statements are about the oracle's integer definitions (`lean/oracle/DdgOracle/NormalCoords/
Formulas.lean`, namespace `DDG.NormalCoords.Formulas`), i.e. about the code the oracle runs.

Face convention: side `k` runs from corner `k` to corner `k + 1`; corner `k` is opposite side
`k + 1`. A face `(n₀, n₁, n₂)` has the corner coordinates
`c_k = cornerCoord n_{k+1} n_k n_{k+2}` and the emanating counts `e_k = emanating n_{k+1} n_k n_{k+2}`.

* **`ic-normal-sign`** — `n⁺ − n⁻ = n`, `n⁺ n⁻ = 0`, both nonnegative, `n⁻ = 1` exactly on
  `n = −1` (`nPlus_sub_nMinus`, `nPlus_mul_nMinus`, `nMinus_eq_one_iff`).
* **`ic-emanating`** (arithmetic half) — `e ≥ 0`, and at most one corner of a face has `e > 0`
  (`emanating_exclusive`): item 2 of SPEC §2.4 is implied by Eq. 2 itself.
* **`ic-corner-coordinates`** (arithmetic half) — the numerator of Eq. 1 is never negative
  (`cornerNumerator_nonneg`), so the `/ 2` of the oracle never sees a negative integer (SPEC I14);
  on counts of the form `n_k = c_k + c_{k+1} + e_{k+2}` with at most one `e > 0` and `c = 0` at an
  emanating corner, Eq. 1 and 2 return `c` and `e` (`corner_eval`).
* **`ic-face-validity`** — the face is valid (the three numerators even) if and only if
  `n⁺₀ + n⁺₁ + n⁺₂` is even or some corner emanates (`faceValid_iff`), and then it has the form above
  (`exists_params_of_faceValid`).
* **`ic-reverse-crossing`** — reversal is an involution on `0 ≤ p < n` (`reverseCrossing_reverseCrossing`,
  `reverseCrossing_mem`).
* **`ic-flip-normal`** (identity half) — on two faces given by their arcs, every corner quantity of
  Eq. 3 is the arc count it names (`flipDoubled_eval`); `DDG/NormalCoords/Flip.lean` turns
  this into the count of arcs crossing the new diagonal and proves **`ic-flip-involution`**.
* **`ic-split-edge-normal`** — Eq. 6 counts every arc of a face with a side carrying no
  crossing (`splitEdgeNormal_eq`).

Corrections to the paper / SPEC found here: Eq. 3's first term is printed `c^{jk}_l`, which names no
corner; it is `c^{ji}_l` (corner `l` of `jil`), which is what `DDG.NormalCoords.Formulas.flipDoubled`
computes. Eq. 6 is `max(n_ki, n_jk, 0)` (SPEC §4.5 drops the `0`, wrong when both are `−1`).
-/

namespace DDG

namespace NormalCoords

open DDG.NormalCoords.Formulas

/-! ## Signs -/

/-- `n⁺ − n⁻ = n`. [formula: ic-normal-sign] -/
theorem nPlus_sub_nMinus (n : ℤ) : nPlus n - nMinus n = n := by
  simp only [nPlus, nMinus]; omega

/-- `n⁺ n⁻ = 0`. [formula: ic-normal-sign] -/
theorem nPlus_mul_nMinus (n : ℤ) : nPlus n * nMinus n = 0 := by
  rcases le_or_gt 0 n with h | h
  · have : nMinus n = 0 := by simp only [nMinus]; omega
    rw [this, mul_zero]
  · have : nPlus n = 0 := by simp only [nPlus]; omega
    rw [this, zero_mul]

/-- `n⁺ ≥ 0` and `n⁻ ≥ 0`. [formula: ic-normal-sign] -/
theorem nPlus_nonneg (n : ℤ) : 0 ≤ nPlus n ∧ 0 ≤ nMinus n := by
  simp only [nPlus, nMinus]; omega

/-- `n⁻ = 1` exactly on a shared edge `n = −1`. [formula: ic-normal-sign] -/
theorem nMinus_eq_one_iff (n : ℤ) : nMinus n = 1 ↔ n = -1 := by
  simp only [nMinus]; omega

/-- [formula: ic-normal-sign] -/
theorem nMinus_eq_zero_iff (n : ℤ) : nMinus n = 0 ↔ 0 ≤ n := by
  simp only [nMinus]; omega

/-! ## Emanating counts and corner coordinates -/

/-- `e ≥ 0`. [formula: ic-emanating] -/
theorem emanating_nonneg (a b c : ℤ) : 0 ≤ emanating a b c := by
  simp only [emanating]; omega

/-- At most one corner of a face emanates: if `e^{ij}_k > 0` then `e^{jk}_i = e^{ki}_j = 0`.
[formula: ic-emanating] -/
theorem emanating_exclusive {a b c : ℤ} (h : 0 < emanating a b c) :
    emanating b c a = 0 ∧ emanating c a b = 0 := by
  simp only [emanating, nPlus] at *; omega

/-- `emanating` is symmetric in the two sides at the corner. [formula: ic-emanating] -/
theorem emanating_comm (a b c : ℤ) : emanating a b c = emanating a c b := by
  simp only [emanating]; omega

/-- The numerator of Eq. 1 is never negative. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_nonneg (a b c : ℤ) : 0 ≤ cornerNumerator a b c := by
  simp only [cornerNumerator, emanating, nPlus]; omega

/-- `cornerNumerator` is symmetric in the two sides at the corner. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_comm (a b c : ℤ) : cornerNumerator a b c = cornerNumerator a c b := by
  simp only [cornerNumerator, emanating, nPlus]; omega

/-- An emanating corner is cut by no arc: `e^{ij}_k > 0 → c^{ij}_k = 0`. [formula: ic-corner-coordinates] -/
theorem cornerNumerator_eq_zero_of_emanating {a b c : ℤ} (h : 0 < emanating a b c) :
    cornerNumerator a b c = 0 := by
  simp only [cornerNumerator, emanating, nPlus] at *; omega

/-- The checked form of Eq. 1 agrees with the plain one exactly when the numerator is even
(it is never negative). [formula: ic-corner-coordinates] -/
theorem cornerCoord?_eq (a b c : ℤ) :
    cornerCoord? a b c = if cornerNumerator a b c % 2 = 0 then some (cornerCoord a b c) else none := by
  have h := cornerNumerator_nonneg a b c
  simp only [cornerCoord?, cornerCoord]
  by_cases h2 : cornerNumerator a b c % 2 = 0
  · simp [h2, not_lt.mpr h]
  · simp [h2, not_lt.mpr h]

/-- **Eq. 1 and 2 on a face given by its arcs.** If the side counts (as `n⁺`) are
`n₀ = c₀ + c₁ + e₂`, `n₁ = c₁ + c₂ + e₀`, `n₂ = c₂ + c₀ + e₁` (side `k` carries the arcs that cut
its two corners and those from the opposite corner) with at most one `e > 0` and `c = 0` at an
emanating corner, then Eq. 1 at corner `0` is `2c₀` before halving and Eq. 2 at corner `0` is `e₀`;
by symmetry the same at every corner. [formula: ic-corner-coordinates] [formula: ic-emanating] -/
theorem corner_eval {n0 n1 n2 c0 c1 c2 e0 e1 e2 : ℤ} (hc0 : 0 ≤ c0) (hc1 : 0 ≤ c1) (hc2 : 0 ≤ c2)
    (he0 : 0 ≤ e0) (he1 : 0 ≤ e1) (he2 : 0 ≤ e2)
    (hx : (e0 = 0 ∧ e1 = 0) ∨ (e1 = 0 ∧ e2 = 0) ∨ (e2 = 0 ∧ e0 = 0))
    (hz0 : e0 = 0 ∨ c0 = 0) (hz1 : e1 = 0 ∨ c1 = 0) (hz2 : e2 = 0 ∨ c2 = 0)
    (h0 : nPlus n0 = c0 + c1 + e2) (h1 : nPlus n1 = c1 + c2 + e0) (h2 : nPlus n2 = c2 + c0 + e1) :
    cornerNumerator n1 n0 n2 = 2 * c0 ∧ emanating n1 n0 n2 = e0 ∧ cornerCoord n1 n0 n2 = c0 := by
  have hn : cornerNumerator n1 n0 n2 = 2 * c0 := by
    simp only [cornerNumerator, emanating]; rw [h0, h1, h2]; omega
  refine ⟨hn, ?_, ?_⟩
  · simp only [emanating]; rw [h0, h1, h2]; omega
  · simp only [cornerCoord]; rw [hn]; omega

/-- **Validity of a face** (SPEC §2.4 item 1): the three numerators of Eq. 1 are even if and only
if `n⁺₀ + n⁺₁ + n⁺₂` is even or one corner emanates. [formula: ic-face-validity] -/
theorem faceValid_iff (n0 n1 n2 : ℤ) :
    faceValid n0 n1 n2 = true ↔
      (nPlus n0 + nPlus n1 + nPlus n2) % 2 = 0 ∨ 0 < emanating n1 n0 n2 ∨
        0 < emanating n2 n1 n0 ∨ 0 < emanating n0 n2 n1 := by
  simp only [faceValid, cornerNumerator, emanating, nPlus, Bool.and_eq_true, beq_iff_eq]
  omega

/-- **A valid face has the form of `corner_eval`**, with `c` and `e` given by Eq. 1 and 2.
[formula: ic-face-validity] [formula: ic-corner-coordinates] -/
theorem exists_params_of_faceValid {n0 n1 n2 : ℤ} (h : faceValid n0 n1 n2 = true) :
    nPlus n0 = cornerCoord n1 n0 n2 + cornerCoord n2 n1 n0 + emanating n0 n2 n1 ∧
      nPlus n1 = cornerCoord n2 n1 n0 + cornerCoord n0 n2 n1 + emanating n1 n0 n2 ∧
      nPlus n2 = cornerCoord n0 n2 n1 + cornerCoord n1 n0 n2 + emanating n2 n1 n0 := by
  simp only [faceValid, cornerCoord, cornerNumerator, emanating, nPlus, Bool.and_eq_true,
    beq_iff_eq] at *
  omega

/-! ## Faces given by their arcs -/

/-- The corner data of a face (`c_A, c_B, c_C`, `e_A, e_B, e_C`) is admissible: nonnegative, at most
one corner emanates, and no arc cuts an emanating corner. Sides: `AB` carries `c_A + c_B + e_C`
crossings, and so on. [formula: ic-corner-coordinates] -/
def FaceParams (cA cB cC eA eB eC : ℤ) : Prop :=
  0 ≤ cA ∧ 0 ≤ cB ∧ 0 ≤ cC ∧ 0 ≤ eA ∧ 0 ≤ eB ∧ 0 ≤ eC ∧
    ((eA = 0 ∧ eB = 0 ∧ (eC = 0 ∨ cC = 0)) ∨ (eB = 0 ∧ eC = 0 ∧ (eA = 0 ∨ cA = 0)) ∨
      (eC = 0 ∧ eA = 0 ∧ (eB = 0 ∨ cB = 0)))

/-- Relabelling the corners keeps admissibility. -/
theorem FaceParams.swap {cA cB cC eA eB eC : ℤ} (h : FaceParams cA cB cC eA eB eC) :
    FaceParams cB cA cC eB eA eC := by unfold FaceParams at *; omega

/-- Relabelling the corners keeps admissibility. -/
theorem FaceParams.rot {cA cB cC eA eB eC : ℤ} (h : FaceParams cA cB cC eA eB eC) :
    FaceParams cB cC cA eB eC eA := by unfold FaceParams at *; omega

/-- Eq. 2 on a face given by its arcs returns the arcs from the corner. [formula: ic-emanating] -/
theorem emanating_eval {o x y cA cB cC eA eB eC : ℤ} (hp : FaceParams cA cB cC eA eB eC)
    (ho : nPlus o = cA + cB + eC) (hx : nPlus x = cC + cA + eB) (hy : nPlus y = cB + cC + eA) :
    emanating o x y = eC := by
  obtain ⟨h1, h2, h3, h4, h5, h6, hc⟩ := hp
  simp only [emanating]
  rw [ho, hx, hy]
  omega

/-- Eq. 1 (before halving) on a face given by its arcs is twice the arcs around the corner.
[formula: ic-corner-coordinates] -/
theorem numerator_eval {o x y cA cB cC eA eB eC : ℤ} (hp : FaceParams cA cB cC eA eB eC)
    (ho : nPlus o = cA + cB + eC) (hx : nPlus x = cC + cA + eB) (hy : nPlus y = cB + cC + eA) :
    cornerNumerator o x y = 2 * cC := by
  obtain ⟨h1, h2, h3, h4, h5, h6, hc⟩ := hp
  simp only [cornerNumerator, emanating]
  rw [ho, hx, hy]
  omega

/-- [formula: ic-corner-coordinates] -/
theorem cornerCoord_of_numerator {o x y cC : ℤ} (hn : cornerNumerator o x y = 2 * cC) :
    cornerCoord o x y = cC := by
  simp only [cornerCoord]; rw [hn]; omega

/-- **Eq. 1 and 2 on a face given by its arcs**, at corner `C` (opposite the side `o` between `A`
and `B`; `x` the side between `C` and `A`, `y` between `B` and `C`).
[formula: ic-corner-coordinates] [formula: ic-emanating] -/
theorem coord_eval {o x y cA cB cC eA eB eC : ℤ} (hp : FaceParams cA cB cC eA eB eC)
    (ho : nPlus o = cA + cB + eC) (hx : nPlus x = cC + cA + eB) (hy : nPlus y = cB + cC + eA) :
    cornerCoord o x y = cC ∧ emanating o x y = eC :=
  ⟨cornerCoord_of_numerator (numerator_eval hp ho hx hy), emanating_eval hp ho hx hy⟩

/-- [formula: ic-corner-coordinates] -/
theorem cornerCoord_nonneg (a b c : ℤ) : 0 ≤ cornerCoord a b c := by
  have := cornerNumerator_nonneg a b c
  simp only [cornerCoord]; omega

/-- **A valid face is a face given by its arcs**, with `c`, `e` from Eq. 1, 2.
[formula: ic-face-validity] -/
theorem params_of_faceValid {n0 n1 n2 : ℤ} (h : faceValid n0 n1 n2 = true) :
    ∃ c0 c1 c2 e0 e1 e2, FaceParams c0 c1 c2 e0 e1 e2 ∧ nPlus n0 = c0 + c1 + e2 ∧
      nPlus n1 = c1 + c2 + e0 ∧ nPlus n2 = c2 + c0 + e1 := by
  refine ⟨cornerCoord n1 n0 n2, cornerCoord n2 n1 n0, cornerCoord n0 n2 n1, emanating n1 n0 n2,
    emanating n2 n1 n0, emanating n0 n2 n1, ?_, exists_params_of_faceValid h⟩
  have k0 : emanating n1 n0 n2 = 0 ∨
      (cornerCoord n1 n0 n2 = 0 ∧ emanating n0 n2 n1 = 0 ∧ emanating n2 n1 n0 = 0) := by
    rcases (emanating_nonneg n1 n0 n2).eq_or_lt with h0 | h0
    · exact Or.inl h0.symm
    · have := emanating_exclusive h0
      refine Or.inr ⟨?_, this.1, this.2⟩
      simp only [cornerCoord, cornerNumerator_eq_zero_of_emanating h0]; rfl
  have k1 : emanating n2 n1 n0 = 0 ∨
      (cornerCoord n2 n1 n0 = 0 ∧ emanating n1 n0 n2 = 0 ∧ emanating n0 n2 n1 = 0) := by
    rcases (emanating_nonneg n2 n1 n0).eq_or_lt with h0 | h0
    · exact Or.inl h0.symm
    · have := emanating_exclusive h0
      refine Or.inr ⟨?_, this.1, this.2⟩
      simp only [cornerCoord, cornerNumerator_eq_zero_of_emanating h0]; rfl
  have k2 : emanating n0 n2 n1 = 0 ∨
      (cornerCoord n0 n2 n1 = 0 ∧ emanating n2 n1 n0 = 0 ∧ emanating n1 n0 n2 = 0) := by
    rcases (emanating_nonneg n0 n2 n1).eq_or_lt with h0 | h0
    · exact Or.inl h0.symm
    · have := emanating_exclusive h0
      refine Or.inr ⟨?_, this.1, this.2⟩
      simp only [cornerCoord, cornerNumerator_eq_zero_of_emanating h0]; rfl
  have g0 := cornerCoord_nonneg n1 n0 n2
  have g1 := cornerCoord_nonneg n2 n1 n0
  have g2 := cornerCoord_nonneg n0 n2 n1
  have f0 := emanating_nonneg n1 n0 n2
  have f1 := emanating_nonneg n2 n1 n0
  have f2 := emanating_nonneg n0 n2 n1
  clear h
  generalize cornerCoord n1 n0 n2 = c0 at *
  generalize cornerCoord n2 n1 n0 = c1 at *
  generalize cornerCoord n0 n2 n1 = c2 at *
  generalize emanating n1 n0 n2 = e0 at *
  generalize emanating n2 n1 n0 = e1 at *
  generalize emanating n0 n2 n1 = e2 at *
  unfold FaceParams
  omega

/-! ## Edge flip -/

/-- **Eq. 3 on two faces given by their arcs**: every corner quantity of Eq. 3 is the arc count it
names. [formula: ic-flip-normal] -/
theorem flipDoubled_eval {nij njk nki nil nlj ci cj ck ei ej ek cj' ci' cl ej' ei' el : ℤ}
    (p1 : FaceParams ci cj ck ei ej ek) (p2 : FaceParams cj' ci' cl ej' ei' el)
    (h_ij : nPlus nij = ci + cj + ek) (h_jk : nPlus njk = cj + ck + ei)
    (h_ki : nPlus nki = ck + ci + ej) (h_ij' : nPlus nij = cj' + ci' + el)
    (h_il : nPlus nil = ci' + cl + ej') (h_lj : nPlus nlj = cl + cj' + ei') :
    flipDoubled nij njk nki nil nlj =
      2 * cl + 2 * ck + iabs (cj' - cj) + iabs (ci' - ci) - el - ek + 2 * (ei' + ei + ej' + ej) := by
  have q1 : FaceParams ci' cj' cl ei' ej' el := p2.swap
  have q2 : FaceParams ci' cl cj' ei' el ej' := p2.rot
  have q3 : FaceParams cj' cl ci' ej' el ei' := p2.rot.rot.swap
  have q4 : FaceParams cj ci ck ej ei ek := p1.swap
  have q5 : FaceParams ci ck cj ei ek ej := p1.rot.rot.swap
  have q6 : FaceParams cj ck ci ej ek ei := p1.rot
  have l1 := coord_eval (o := nij) (x := nil) (y := nlj) q1 (by linarith) (by linarith) (by linarith)
  have l2 := coord_eval (o := nil) (x := nij) (y := nlj) q2 (by linarith) (by linarith) (by linarith)
  have l3 := coord_eval (o := nlj) (x := nij) (y := nil) q3 (by linarith) (by linarith) (by linarith)
  have r1 := coord_eval (o := nij) (x := njk) (y := nki) q4 (by linarith) (by linarith) (by linarith)
  have r2 := coord_eval (o := nki) (x := nij) (y := njk) q5 (by linarith) (by linarith) (by linarith)
  have r3 := coord_eval (o := njk) (x := nij) (y := nki) q6 (by linarith) (by linarith) (by linarith)
  simp only [flipDoubled]
  rw [l1.1, l2.1, l3.1, r1.1, r2.1, r3.1, l1.2, l2.2, l3.2, r1.2, r2.2, r3.2]

/-! ## Crossings -/

/-- Reversal of a combinatorial crossing is an involution. [formula: ic-reverse-crossing] -/
theorem reverseCrossing_reverseCrossing (n p : ℤ) : reverseCrossing n (reverseCrossing n p) = p := by
  simp only [reverseCrossing]; omega

/-- Reversal maps the crossings `0 ≤ p < n` of a halfedge onto those of its twin, reversing the
order. [formula: ic-reverse-crossing] -/
theorem reverseCrossing_mem {n p : ℤ} : (0 ≤ reverseCrossing n p ∧ reverseCrossing n p < n) ↔
    (0 ≤ p ∧ p < n) := by
  simp only [reverseCrossing]; omega

/-- [formula: ic-reverse-crossing] -/
theorem reverseCrossing_lt_iff (n p q : ℤ) : reverseCrossing n p < reverseCrossing n q ↔ q < p := by
  simp only [reverseCrossing]; omega

/-! ## Edge split along a shared edge -/

/-- **Eq. 6** counts every arc of the face `ijk` when the side `ij` carries no crossing: with
`n⁺_ij = 0` the face has `c_i = c_j = e_k = 0`, its arcs are the `c_k` arcs around `k` and the
`e_i + e_j` arcs from `i` or `j` (at most one of them nonzero), and `max(n_ki, n_jk, 0)` is their
number. Arguments: `n_ij`, `n_jk`, `n_ki` as `n₀, n₁, n₂`. [formula: ic-split-edge-normal] -/
theorem splitEdgeNormal_eq {n0 n1 n2 : ℤ} (hv : faceValid n0 n1 n2 = true) (h0 : nPlus n0 = 0) :
    splitEdgeNormal n2 n1 = cornerCoord n0 n2 n1 + emanating n1 n0 n2 + emanating n2 n1 n0 ∧
      cornerCoord n1 n0 n2 = 0 ∧ cornerCoord n2 n1 n0 = 0 ∧ emanating n0 n2 n1 = 0 := by
  simp only [faceValid, splitEdgeNormal, cornerCoord, cornerNumerator, emanating, nPlus,
    Bool.and_eq_true, beq_iff_eq] at *
  omega

end NormalCoords

end DDG
