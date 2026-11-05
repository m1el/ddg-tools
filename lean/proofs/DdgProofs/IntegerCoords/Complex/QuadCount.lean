import DdgProofs.IntegerCoords.Complex.Quad

/-!
# `ic-flip-normal` (the C part): Eq. 3 after gluing, and the arcs at the four corners

(`Hyp.outer`, the curves, is in `Complex/QuadFlip.lean`.)

`Complex/Quad.lean` explains the setting. Here:

* `Hyp.flipNormal_eq`: Eq. 3 as the oracle computes it (`flipNormal`, the corrected first terms),
  from the five signed normal coordinates of the two faces, is `d − D_kl`: the number of arcs of
  the glued pair crossing `kl`, minus the arcs from `k` to `l`;
* `Hyp.outer`: the re-cut arcs are the same curves, on all four outer sides;
* `Hyp.em_k`, `Hyp.em_l`, `Hyp.em_i`, `Hyp.em_j`: the arcs leaving the corners. The old corner `k`
  of `ijk` is cut by `kl` into the corner `k` of `F₁ = lki` and the corner `k` of `F₂ = klj`, with
  `n⁻_kl = D_kl` between them; the new corner `i` of `F₁` holds the arcs of the corners `i` of
  `jil` and of `ijk` and `n⁻_ij` between them (and the same at `l`, `j`).
-/

namespace DDG.IntegerCoords

open Ddg.IntegerCoords (cornerCoord emanating faceValid nPlus nMinus flipNormal flipDoubled iabs)

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

/-- The relative position of the blocks `[x0, x0 + a2)` (arcs from `k`) and `[y1, y1 + b2)` (arcs
from `l`) on the diagonal, which makes `P1`, `P2`, `D_kl` linear. -/
lemma kinds_cases {x0 a2 y1 b2 p1 p2 D : ℕ} (dP1 : p1 = y1 - (x0 + a2)) (dP2 : p2 = x0 - (y1 + b2))
    (dD : D = min (x0 + a2) (y1 + b2) - max x0 y1) :
    (x0 + a2 ≤ y1 ∧ p1 + x0 + a2 = y1 ∧ p2 = 0 ∧ D = 0) ∨
    (y1 + b2 ≤ x0 ∧ p2 + y1 + b2 = x0 ∧ p1 = 0 ∧ D = 0) ∨
    (y1 < x0 + a2 ∧ x0 < y1 + b2 ∧ p1 = 0 ∧ p2 = 0 ∧ D + max x0 y1 = min (x0 + a2) (y1 + b2)) := by
  omega

set_option hygiene false in
/-- `glue_facts`, then `P1`, `P2`, `D_kl` made linear by `kinds_cases`. -/
macro "glue_lin" h:ident : tactic =>
  `(tactic| (glue_facts $h <;> rcases kinds_cases dP1 dP2 dD with ⟨k1, k2, k3, k4⟩ | ⟨k1, k2, k3, k4⟩ | ⟨k1, k2, k3, k4, k5⟩ <;> clear dP1 dP2 dD))

lemma x_eq (k : Fin 3) : cc nA k = x nA k := rfl
lemma y_eq (k : Fin 3) : cc nB k = y nB k := rfl
lemma ea_eq (k : Fin 3) : em nA (cc nA) k = ea nA k := rfl
lemma eb_eq (k : Fin 3) : em nB (cc nB) k = eb nB k := rfl

set_option maxHeartbeats 1000000 in
/-- **Eq. 3 counts the arcs of the glued pair crossing `kl`** (`ic-flip-normal`, the C part). For two
valid faces `ijk`, `jil` on `ij` (`Hyp`), the oracle's `flipNormal` of the signed normal
coordinates `n_ij = n⁺_ij − n⁻_ij`, `n_jk`, `n_ki`, `n_il`, `n_lj` (the others enter only through
`n⁺`) is `d − D_kl`: `d` arcs cross the new diagonal, and `D_kl` arcs run from `k` to `l` (with
one, the new edge runs along that input edge: `n_kl = −1`; with two or more, which no straight
arcs can do, the code refuses `n_kl < −1`). [formula: ic-flip-normal] -/
theorem Hyp.flipNormal_eq (h : Hyp nA nB m) {nij njk nki nil nlj : ℤ}
    (hij : nij = (nA 0 : ℤ) - m) (hjk : nPlus njk = nA 1) (hki : nPlus nki = nA 2)
    (hil : nPlus nil = nB 1) (hlj : nPlus nlj = nB 2) :
    flipNormal nij njk nki nil nlj = (d nA nB m : ℤ) - Dkl nA nB := by
  have hA := admissible_cc h.vA
  have hB := admissible_cc h.vB
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  have p1 := hA.faceParams 0
  have p2 := hB.faceParams 0
  simp only [l1, l2, x_eq (nA := nA), y_eq (nB := nB), ea_eq (nA := nA), eb_eq (nB := nB)] at p1 p2
  have hm := h.hm
  have hmN := h.hmN
  have hpij : nPlus nij = nA 0 := by
    simp only [nPlus, hij]; rcases Nat.lt_or_ge 0 m with h0 | h0
    · have := hmN (by omega); omega
    · omega
  have hmij : nMinus nij = m := by
    simp only [nMinus, hij]; rcases Nat.lt_or_ge 0 m with h0 | h0
    · have := hmN (by omega); omega
    · omega
  obtain ⟨fA0, fA1, fA2, -, -, -, fB0, fB1, fB2, -, -, -⟩ := h.facts
  have e := flipDoubled_eval (nij := nij) (njk := njk) (nki := nki) (nil := nil) (nlj := nlj) p1 p2
    (by rw [hpij]; exact_mod_cast fA0) (by rw [hjk]; exact_mod_cast fA1)
    (by rw [hki]; push_cast; omega) (by rw [hpij, h.hN]; exact_mod_cast fB0)
    (by rw [hil]; exact_mod_cast fB1) (by rw [hlj]; push_cast; omega)
  rw [flipNormal, e, hmij]
  simp only [iabs]
  glue_lin h
  all_goals (split_ifs <;> push_cast <;> omega)

/-- The corner `k` of `ijk` is split by `kl`: its arcs are those of corner `k` of `F₁ = lki` (to
`il`), the `D_kl` arcs to `l` (the new edge when it runs along an input edge), and those of corner
`k` of `F₂ = klj` (to `lj`), in this counter-clockwise order. -/
theorem Hyp.em_k (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    ea nA 2 = em (n1 nA nB m) (c1 nA nB) 1 + Dkl nA nB + em (n2 nA nB m) (c2 nA nB) 0 := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [em, n1, n2, c1, c2, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  glue_lin h
  all_goals omega

/-- The corner `l` of `jil` is split by `kl` (corner `l` of `F₁`, the `D_kl` arcs, corner `l` of
`F₂`). -/
theorem Hyp.em_l (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    eb nB 2 = em (n1 nA nB m) (c1 nA nB) 0 + Dkl nA nB + em (n2 nA nB m) (c2 nA nB) 1 := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [em, n1, n2, c1, c2, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  glue_lin h
  all_goals omega

/-- The new corner `i` of `F₁` holds the arcs of corner `i` of `jil`, the old edge `ij` when it runs
along an input edge (`m = n⁻_ij`), and the arcs of corner `i` of `ijk`. -/
theorem Hyp.em_i (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    em (n1 nA nB m) (c1 nA nB) 2 = eb nB 1 + m + ea nA 0 := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [em, n1, c1, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  glue_lin h
  all_goals omega

/-- The new corner `j` of `F₂` holds the arcs of corner `j` of `ijk`, `n⁻_ij`, and those of corner `j`
of `jil`. -/
theorem Hyp.em_j (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    em (n2 nA nB m) (c2 nA nB) 2 = ea nA 1 + m + eb nB 0 := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  simp only [em, n2, c2, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  glue_lin h
  all_goals omega

end Glue

end DDG.IntegerCoords
