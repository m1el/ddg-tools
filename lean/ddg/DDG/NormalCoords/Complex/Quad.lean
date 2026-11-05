import DDG.NormalCoords.Complex.Basic
import DDG.NormalCoords.Formulas
import Mathlib.Data.Fin.VecNotation

/-!
# Integer coordinates, tier C: gluing two face models and re-cutting (`ic-flip-normal`)

`crates/ddg-integer-coords/SPEC.md` §4.2, §9.2 (`ic-flip-normal`, the C part); the paper Sec. 3.2,
Eq. 3; the Rust `flipped_integers` (`crates/ddg-integer-coords/src/flip.rs`).

The edge `ij` with the faces `ijk` and `jil`, in the face-model numbering with the edge on side `0`
of both: face `A = ijk` (corners `0, 1, 2` = `i, j, k`; side `0` = `i→j`, `1` = `j→k`, `2` = `k→i`)
and face `B = jil` (corners `j, i, l`; side `0` = `j→i`, `1` = `i→l`, `2` = `l→j`), with the same
`n⁺_ij` on side `0`, and `m = n⁻_ij` (`m = 1`: `ij` runs along an input edge, then `n⁺_ij = 0`).
Crossing `p` of `A`'s side `0` is crossing `n⁺_ij − 1 − p` of `B`'s (`ic-reverse-crossing`). The
quadrilateral `i, l, j, k` has the four **outer sides** `jk`, `ki` (of `A`) and `il`, `lj` (of `B`),
each with its own direction and numbering, which the flip does not change (the Rust flip keeps
every other edge's id, direction and normal coordinate).

After the flip the faces are `F₁ = lki` (side `0` = `l→k`, `1` = `k→i`, `2` = `i→l`) and
`F₂ = klj` (side `0` = `k→l`, `1` = `l→j`, `2` = `j→k`) — the Rust faces `(i, l, k)` and `(j, k, l)`
with the new edge moved to side `0` — with the side counts `(d, n_ki, n_il)` and `(d, n_lj, n_jk)`,
`d = n⁺_kl`.

* `Hyp.flipNormal_eq` (`QuadCount.lean`): Eq. 3 (the oracle's `flipNormal`, with the corrected
  first terms) of the two faces is `d − D_kl`, where `D_kl` is the number of arcs from `k` to `l`
  through `ij` (`D_kl ≤ 1` ⟺ `n_kl ≥ −1`, the Rust refusal of `n_kl < −1`) and `d` the explicit
  count `Glue.d` of the arcs crossing `kl`.
* `Hyp.admissible1/2`, `Hyp.cc1/2` (`QuadAdm1.lean`, `QuadAdm2.lean`): the new faces are valid, with
  the corner counts of the re-cut (`c1`, `c2`) as their Eq. 1, so the new face models are the
  canonical ones of their counts.
* **`Hyp.outer`** (`QuadFlip.lean`; one file per side, `QuadOuterJK.lean` …, so that the case
  analyses build in parallel), the correctness of Eq. 3 after gluing: for every crossing of every
  outer side, the curve that enters the quadrilateral there leaves it at the same place (crossing
  of an outer side, or a corner `i, j, k, l`) whether it is followed through the old faces (across
  `ij`) or through the new ones (across `kl`). So the re-cut arcs are the same curves.
* `Hyp.em_k`, `Hyp.em_l`, `Hyp.em_i`, `Hyp.em_j` (`QuadCount.lean`): the arcs leaving the corners:
  the old corner `k` (`l`) is split by `kl` into two corners whose arcs, with the `D_kl` arcs to `l`
  between them (`n⁻_kl`), are the old ones; the new corner `i` (`j`) holds the arcs of the two old
  corners and `n⁻_ij` between them (used for Eq. 4, `Complex/Flip.lean`).

**Build time.** `omega` is fast on these goals only when the facts are linear and contain no
implication: `glue_facts` names the atoms, splits by the emanating corner of each face (16 cases)
and the definitions of `D_kl` (`min`/`max`) are cleared or made linear before `omega`
(`QuadCount.kinds_cases`).
-/

namespace DDG.IntegerCoords

open DDG.IntegerCoords.Formulas (cornerCoord emanating faceValid nPlus nMinus flipNormal flipDoubled iabs)

namespace Glue

variable (nA nB : Fin 3 → ℕ) (m : ℕ)

/-- Corner counts (Eq. 1) of `A = ijk` (`x 0 = c_i`, `x 1 = c_j`, `x 2 = c_k`). -/
def x (k : Fin 3) : ℕ := cc nA k
/-- Corner counts (Eq. 1) of `B = jil` (`y 0 = c_j`, `y 1 = c_i`, `y 2 = c_l`). -/
def y (k : Fin 3) : ℕ := cc nB k
/-- Arcs leaving the corners of `A` (`ea 2 = e_k`, to `ij`). -/
def ea (k : Fin 3) : ℕ := em nA (cc nA) k
/-- Arcs leaving the corners of `B` (`eb 2 = e_l`, to `ij`). -/
def eb (k : Fin 3) : ℕ := em nB (cc nB) k

/-- Arcs from `il` to `jk` (cut `j` in `A`, cut `i` in `B`). -/
def P1 : ℕ := y nB 1 - (x nA 0 + ea nA 2)
/-- Arcs from `lj` to `ki`. -/
def P2 : ℕ := x nA 0 - (y nB 1 + eb nB 2)
/-- Arcs from `k` to `l` through `ij`. -/
def Dkl : ℕ := min (x nA 0 + ea nA 2) (y nB 1 + eb nB 2) - max (x nA 0) (y nB 1)
/-- `n⁺_kl`: the arcs crossing `kl` (corner arcs at `k` and `l`, the arcs between opposite sides,
the arcs from `i` and `j`, and the input edge along `ij`). -/
def d : ℕ := x nA 2 + y nB 2 + P1 nA nB + P2 nA nB + eb nB 1 + ea nA 0 + eb nB 0 + ea nA 1 + m

/-- Side counts of `F₁ = lki`. -/
def n1 : Fin 3 → ℕ := ![d nA nB m, nA 2, nB 1]
/-- Side counts of `F₂ = klj`. -/
def n2 : Fin 3 → ℕ := ![d nA nB m, nB 2, nA 1]
/-- Corner counts of `F₁` (`l`, `k`, `i`). -/
def c1 : Fin 3 → ℕ := ![y nB 2 + eb nB 0 + P1 nA nB, x nA 2 + ea nA 1 + P2 nA nB, min (x nA 0) (y nB 1)]
/-- Corner counts of `F₂` (`k`, `l`, `j`). -/
def c2 : Fin 3 → ℕ := ![x nA 2 + ea nA 0 + P1 nA nB, y nB 2 + eb nB 1 + P2 nA nB, min (x nA 1) (y nB 0)]

/-- The hypotheses: two valid faces on `ij` with the same `n⁺_ij`, `m = n⁻_ij`. -/
structure Hyp : Prop where
  vA : faceValid (nA 0) (nA 1) (nA 2) = true
  vB : faceValid (nB 0) (nB 1) (nB 2) = true
  hN : nA 0 = nB 0
  hm : m ≤ 1
  hmN : m = 1 → nA 0 = 0

variable {nA nB m}

/-- The relations between side counts, corner counts and emanating counts of the two faces. -/
lemma Hyp.facts (h : Hyp nA nB m) :
    nA 0 = x nA 0 + x nA 1 + ea nA 2 ∧ nA 1 = x nA 1 + x nA 2 + ea nA 0 ∧
    nA 2 = x nA 2 + x nA 0 + ea nA 1 ∧
    (0 < ea nA 0 → x nA 0 = 0 ∧ ea nA 1 = 0 ∧ ea nA 2 = 0) ∧
    (0 < ea nA 1 → x nA 1 = 0 ∧ ea nA 2 = 0 ∧ ea nA 0 = 0) ∧
    (0 < ea nA 2 → x nA 2 = 0 ∧ ea nA 0 = 0 ∧ ea nA 1 = 0) ∧
    nB 0 = y nB 0 + y nB 1 + eb nB 2 ∧ nB 1 = y nB 1 + y nB 2 + eb nB 0 ∧
    nB 2 = y nB 2 + y nB 0 + eb nB 1 ∧
    (0 < eb nB 0 → y nB 0 = 0 ∧ eb nB 1 = 0 ∧ eb nB 2 = 0) ∧
    (0 < eb nB 1 → y nB 1 = 0 ∧ eb nB 2 = 0 ∧ eb nB 0 = 0) ∧
    (0 < eb nB 2 → y nB 2 = 0 ∧ eb nB 0 = 0 ∧ eb nB 1 = 0) := by
  have hA := admissible_cc h.vA
  have hB := admissible_cc h.vB
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have a0 := hA.n_eq 0; have a1 := hA.n_eq 1; have a2 := hA.n_eq 2
  have b0 := hB.n_eq 0; have b1 := hB.n_eq 1; have b2 := hB.n_eq 2
  have p0 := hA.2.1 0; have p1 := hA.2.1 1; have p2 := hA.2.1 2
  have q0 := hA.2.2 0; have q1 := hA.2.2 1; have q2 := hA.2.2 2
  have r0 := hB.2.1 0; have r1 := hB.2.1 1; have r2 := hB.2.1 2
  have s0 := hB.2.2 0; have s1 := hB.2.2 1; have s2 := hB.2.2 2
  simp only [l1, l2, l3, l4, l5, l6] at *
  simp only [x, y, ea, eb]
  refine ⟨a0, a1, a2, fun h => ⟨p0 h, q0 h⟩, fun h => ⟨p1 h, q1 h⟩, fun h => ⟨p2 h, q2 h⟩,
    b0, b1, b2, fun h => ⟨r0 h, s0 h⟩, fun h => ⟨r1 h, s1 h⟩, fun h => ⟨r2 h, s2 h⟩⟩

/-- At most one corner of a face emanates, and it is cut by no arc. -/
lemma em_cases {x0 x1 x2 a0 a1 a2 : ℕ} (e0 : 0 < a0 → x0 = 0 ∧ a1 = 0 ∧ a2 = 0)
    (e1 : 0 < a1 → x1 = 0 ∧ a2 = 0 ∧ a0 = 0) (e2 : 0 < a2 → x2 = 0 ∧ a0 = 0 ∧ a1 = 0) :
    (a0 = 0 ∧ a1 = 0 ∧ a2 = 0) ∨ (x0 = 0 ∧ a1 = 0 ∧ a2 = 0) ∨ (x1 = 0 ∧ a0 = 0 ∧ a2 = 0) ∨
      (x2 = 0 ∧ a0 = 0 ∧ a1 = 0) := by
  omega

set_option hygiene false in
/-- Introduce the facts of `Hyp.facts` and the definitions of the kinds, as hypotheses about
opaque natural numbers, for `omega` (not hygienic: it names the variables `nA nB m` of the
goal). -/
macro "glue_facts0" h:ident : tactic =>
  `(tactic| (
    obtain ⟨fA0, fA1, fA2, eA0, eA1, eA2, fB0, fB1, fB2, eB0, eB1, eB2⟩ := ($h).facts
    have fN := ($h).hN
    have fm := ($h).hm
    have fmN := ($h).hmN
    have dP1 : P1 nA nB = y nB 1 - (x nA 0 + ea nA 2) := rfl
    have dP2 : P2 nA nB = x nA 0 - (y nB 1 + eb nB 2) := rfl
    have dD : Dkl nA nB = min (x nA 0 + ea nA 2) (y nB 1 + eb nB 2) - max (x nA 0) (y nB 1) := rfl
    have dd : d nA nB m = x nA 2 + y nB 2 + P1 nA nB + P2 nA nB + eb nB 1 + ea nA 0 + eb nB 0 +
      ea nA 1 + m := rfl
    generalize x nA 0 = x0 at *
    generalize x nA 1 = x1 at *
    generalize x nA 2 = x2 at *
    generalize y nB 0 = y0 at *
    generalize y nB 1 = y1 at *
    generalize y nB 2 = y2 at *
    generalize ea nA 0 = a0 at *
    generalize ea nA 1 = a1 at *
    generalize ea nA 2 = a2 at *
    generalize eb nB 0 = b0 at *
    generalize eb nB 1 = b1 at *
    generalize eb nB 2 = b2 at *
    generalize P1 nA nB = p1 at *
    generalize P2 nA nB = p2 at *
    generalize Dkl nA nB = D at *
    generalize d nA nB m = dk at *))

set_option hygiene false in
/-- Split the facts of `glue_facts0` by the emanating corner of each face (no implications left). -/
macro "glue_split" : tactic =>
  `(tactic| (
    rcases em_cases eA0 eA1 eA2 with ⟨hA1, hA2, hA3⟩ | ⟨hA1, hA2, hA3⟩ | ⟨hA1, hA2, hA3⟩ |
      ⟨hA1, hA2, hA3⟩ <;>
    rcases em_cases eB0 eB1 eB2 with ⟨hB1, hB2, hB3⟩ | ⟨hB1, hB2, hB3⟩ | ⟨hB1, hB2, hB3⟩ |
      ⟨hB1, hB2, hB3⟩ <;>
    subst hA1 hA2 hA3 hB1 hB2 hB3 <;>
    clear eA0 eA1 eA2 eB0 eB1 eB2))

/-- `glue_facts0`, then `glue_split`. -/
macro "glue_facts" h:ident : tactic => `(tactic| (glue_facts0 $h; glue_split))

/-- With an arc from `k` to `l`, nothing else crosses `kl`: `d = 0`. -/
lemma Hyp.d_eq_zero (h : Hyp nA nB m) (hD : 0 < Dkl nA nB) : d nA nB m = 0 := by
  glue_facts h
  all_goals omega

/-! ## The arcs through the quadrilateral, before and after -/

/-- The outer sides of the quadrilateral, with the directions of the halfedges `j→k`, `k→i`
(of `A`) and `i→l`, `l→j` (of `B`), which the flip keeps. -/
inductive QS where
  | jk | ki | il | lj
  deriving DecidableEq

/-- A boundary point of the quadrilateral: one of its corners, or a crossing of an outer side
(counted from that side's tail). -/
inductive QP where
  | ci | cj | ck | cl
  | side (s : QS) (p : ℕ)
  deriving DecidableEq

/-- A boundary point of `A = ijk` other than its side `ij`. -/
def ofA : Pt → QP
  | .corner k => if k = 0 then .ci else if k = 1 then .cj else .ck
  | .side k p => if k = 1 then .side .jk p else .side .ki p
/-- A boundary point of `B = jil` other than its side `ji`. -/
def ofB : Pt → QP
  | .corner k => if k = 0 then .cj else if k = 1 then .ci else .cl
  | .side k p => if k = 1 then .side .il p else .side .lj p
/-- A boundary point of `F₁ = lki` other than `lk`. -/
def of1 : Pt → QP
  | .corner k => if k = 0 then .cl else if k = 1 then .ck else .ci
  | .side k p => if k = 1 then .side .ki p else .side .il p
/-- A boundary point of `F₂ = klj` other than `kl`. -/
def of2 : Pt → QP
  | .corner k => if k = 0 then .ck else if k = 1 then .cl else .cj
  | .side k p => if k = 1 then .side .lj p else .side .jk p

variable (nA nB m)

/-- Where an arc of `A` goes, continued across `ij` into `B` (crossing `q` of `A`'s side `0` is
crossing `n⁺_ij − 1 − q` of `B`'s). -/
def thruA : Pt → QP
  | .side k q => if k = 0 then ofB ((model nB).partner 0 (nA 0 - 1 - q)) else ofA (.side k q)
  | .corner k => ofA (.corner k)
/-- Where an arc of `B` goes, continued across `ij` into `A`. -/
def thruB : Pt → QP
  | .side k q => if k = 0 then ofA ((model nA).partner 0 (nB 0 - 1 - q)) else ofB (.side k q)
  | .corner k => ofB (.corner k)
/-- Where an arc of `F₁` goes, continued across `kl` into `F₂`. -/
def thru1 : Pt → QP
  | .side k q => if k = 0 then of2 ((model (n2 nA nB m)).partner 0 (d nA nB m - 1 - q))
      else of1 (.side k q)
  | .corner k => of1 (.corner k)
/-- Where an arc of `F₂` goes, continued across `kl` into `F₁`. -/
def thru2 : Pt → QP
  | .side k q => if k = 0 then of1 ((model (n1 nA nB m)).partner 0 (d nA nB m - 1 - q))
      else of2 (.side k q)
  | .corner k => of2 (.corner k)

/-- The curve entering the quadrilateral at crossing `p` of an outer side, followed through the
old faces (across `ij`): where it leaves the quadrilateral. -/
def oldExit : QS → ℕ → QP
  | .jk, p => thruA nA nB ((model nA).partner 1 p)
  | .ki, p => thruA nA nB ((model nA).partner 2 p)
  | .il, p => thruB nA nB ((model nB).partner 1 p)
  | .lj, p => thruB nA nB ((model nB).partner 2 p)

/-- The same through the new faces (across `kl`). -/
def newExit : QS → ℕ → QP
  | .jk, p => thru2 nA nB m ((model (n2 nA nB m)).partner 2 p)
  | .ki, p => thru1 nA nB m ((model (n1 nA nB m)).partner 1 p)
  | .il, p => thru1 nA nB m ((model (n1 nA nB m)).partner 2 p)
  | .lj, p => thru2 nA nB m ((model (n2 nA nB m)).partner 1 p)

/-- The number of crossings of an outer side. -/
def qcount : QS → ℕ
  | .jk => nA 1
  | .ki => nA 2
  | .il => nB 1
  | .lj => nB 2

variable {nA nB m}

lemma modelA : (model nA).partner = canon nA (x nA) := rfl
lemma modelB : (model nB).partner = canon nB (y nB) := rfl

end Glue

end DDG.IntegerCoords
