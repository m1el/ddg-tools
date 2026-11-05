import DdgProofs.IntegerCoords.Complex.Quad

/-! # The new face `F₂` of the flip is valid (`Complex/Quad.lean`); a file of its own so that the
two case analyses build in parallel. -/

namespace DDG.IntegerCoords

open Ddg.IntegerCoords (faceValid)

namespace Glue

variable {nA nB : Fin 3 → ℕ} {m : ℕ}

set_option maxHeartbeats 4000000 in
theorem Hyp.admissible2 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) : Admissible (n2 nA nB m) (c2 nA nB) := by
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
  refine ⟨fun k => ?_, fun k => ?_, fun k => ?_⟩ <;>
    rcases fin3_cases k with rfl | rfl | rfl <;>
    simp only [em, n2, c2, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
    glue_facts h <;> (clear dD; omega)

lemma Hyp.cc2 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) : cc (n2 nA nB m) = c2 nA nB := by
  have ha := h.admissible2 hD
  have hv : faceValid (n2 nA nB m 0) (n2 nA nB m 1) (n2 nA nB m 2) = true :=
    (FaceModel.exists_valid_iff _).1 ⟨_, rfl, canon_valid ha⟩
  funext k
  have e1 := ha.c_eq k
  have e2 := (admissible_cc hv).c_eq k
  exact_mod_cast e2.trans e1.symm

lemma Hyp.model2 (h : Hyp nA nB m) (hD : Dkl nA nB ≤ 1) :
    (model (n2 nA nB m)).partner = canon (n2 nA nB m) (c2 nA nB) := by
  simp only [model, h.cc2 hD]

end Glue

end DDG.IntegerCoords
