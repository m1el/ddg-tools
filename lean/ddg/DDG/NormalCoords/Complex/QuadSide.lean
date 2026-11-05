import DDG.NormalCoords.Complex.QuadAdm1
import DDG.NormalCoords.Complex.QuadAdm2

/-!
# The tactic for one outer side of the glued quadrilateral

`Complex/QuadOuterJK.lean` … `QuadOuterLJ.lean` (one side per file, so that the four case
analyses build in parallel).
-/

namespace DDG.NormalCoords

namespace Glue

set_option hygiene false in
/-- Unfold the two composites at one outer side, split every case (the corner that emanates in
each face, the blocks of the four canonical models) and close each case by arithmetic. Not
hygienic: it uses the hypotheses `h`, `hD` and the variables `nA nB m` of the goal. -/
macro "glue_side" : tactic =>
  `(tactic| (
    obtain ⟨l1, l2, l3, l4, l5, l6⟩ := fin3_lits
    simp only [oldExit, newExit, thruA, thruB, thru1, thru2]
    simp only [h.model1 hD, h.model2 hD, modelA (nA := nA), modelB (nB := nB)]
    simp only [canon, n1, n2, c1, c2, l1, l2, l3, l4, l5, l6, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
    glue_facts h
    all_goals (
      simp only [Nat.add_zero, Nat.zero_add] at *
      split_ifs <;>
      simp only [ofA, ofB, of1, of2, l1, l2, l3, l4, l5, l6, ite_true, ite_false, Fin.reduceEq,
        QP.side.injEq, reduceCtorEq, true_and, and_true, false_and, and_false] <;>
      try split_ifs <;>
      simp only [ofA, ofB, of1, of2, l1, l2, l3, l4, l5, l6, ite_true, ite_false, Fin.reduceEq,
        QP.side.injEq, reduceCtorEq, true_and, and_true, false_and, and_false] <;>
      try split_ifs <;>
      simp only [ofA, ofB, of1, of2, l1, l2, l3, l4, l5, l6, ite_true, ite_false, Fin.reduceEq,
        QP.side.injEq, reduceCtorEq, true_and, and_true, false_and, and_false] at *)
    all_goals (try simp (disch := omega) only [Nat.sub_sub_self] at *)
    all_goals first | (clear dD; omega) | omega))

end Glue

end DDG.NormalCoords
