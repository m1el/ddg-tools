import DDG.NormalCoords.Complex.Labels
import DDG.NormalCoords.Complex.QuadCount
import DDG.NormalCoords.Complex.QuadAdm1
import DDG.NormalCoords.Complex.QuadAdm2
import Mathlib.Tactic.Ring
import Mathlib.Algebra.Group.Fin.Basic

/-!
# Integer coordinates, tier C: the edge flip of a complex model (Eq. 3, Eq. 4)

`crates/ddg-integer-coords/SPEC.md` §4.2, §9.2 (`ic-flip-normal` C part, `ic-roundabout-update`);
the Rust `flip_edge_with`, `flipped_integers`, `updated_roundabout`
(`crates/ddg-integer-coords/src/flip.rs`) and `Complex::flip_edge` (`ddg-intrinsic`).

**The flip** of the interior edge with sides `(f, s)` (`i→j`, in `ijk`) and `(g, t)` (`j→i`, in
`jil`), `f ≠ g` (an edge whose two sides are in one face has an endpoint of degree one, which
`plan_flip` refuses), exactly as `Complex::flip_edge` renumbers: face `f` becomes `(i, l, k)` with
`i→l` in slot `s`, the new `l→k` in slot `s + 1` and `k→i` kept in slot `s + 2`; face `g` becomes
`(j, k, l)` with `j→k` in slot `t`, `k→l` in slot `t + 1`, `l→j` kept in slot `t + 2`. The
halfedges `j→k` (old `(f, s+1)`) and `i→l` (old `(g, t+1)`) move with their data (`σ`,
`ComplexFlip::move_halfedge_data`); every other side keeps its slot. The new edge has
`n_kl = flipNormal n_ij n_jk n_ki n_il n_lj` (Eq. 3 with the corrected first terms, the Rust
`formulas::flip_normal`).

* `Cx.flip_valid`: if `n_kl ≥ −1` (the Rust refusal of `n_kl < −1`), the flipped complex model is
  valid; `n_kl` is the number of arcs of the glued pair crossing `kl` minus the arcs from `k` to
  `l` (`Glue.Hyp.flipNormal_eq`), and the curves through the quadrilateral are unchanged
  (`Glue.Hyp.outer`).
* `Cx.flip_roundabouts` (**`ic-roundabout-update`**): if the roundabouts of the state are
  consistent (`RoundaboutsOK`, the check of `certify`), so are those after the flip, with `r_kl` and
  `r_lk` by Eq. 4 (the oracle's `roundaboutUpdate` of `r_ki`, the emanating count of corner `k` of
  `kil` and `n⁻_ki`, as `updated_roundabout` computes it; the same at `l`) and every other halfedge
  keeping its roundabout. So Eq. 4 gives the local index of the first input halfedge at or after
  `k→l`: the labels of all ends stay consecutive (`RoundaboutsOK.lab_next`).
-/

set_option linter.unusedSectionVars false

namespace DDG.NormalCoords

open DDG.NormalCoords.Formulas (cornerCoord emanating faceValid nPlus nMinus flipNormal roundaboutUpdate)

/-! ## Rotating a face -/

lemma faceValid_rot (a b c : ℤ) : faceValid a b c = faceValid b c a := by
  simp only [faceValid]
  generalize (DDG.NormalCoords.Formulas.cornerNumerator b a c % 2 == 0) = x
  generalize (DDG.NormalCoords.Formulas.cornerNumerator c b a % 2 == 0) = y
  generalize (DDG.NormalCoords.Formulas.cornerNumerator a c b % 2 == 0) = z
  cases x <;> cases y <;> cases z <;> rfl

lemma cc_rot (n : Fin 3 → ℕ) (r k : Fin 3) : cc (fun k => n (k + r)) k = cc n (k + r) := by
  simp only [cc]
  rw [add_right_comm k 1 r, add_right_comm k 2 r]

lemma em_rot (n : Fin 3 → ℕ) (r k : Fin 3) :
    em (fun k => n (k + r)) (cc fun k => n (k + r)) k = em n (cc n) (k + r) := by
  simp only [em, cc_rot]
  rw [add_right_comm k 1 r, add_right_comm k 2 r]

namespace Cx

variable {F : Type} [DecidableEq F] {C : Cx F}

/-- A flip: the two sides `(f, s)` (`i→j`) and `(g, t)` (`j→i`) of an interior edge, in different
faces. -/
structure FlipAt (C : Cx F) where
  f : F
  g : F
  s : Fin 3
  t : Fin 3
  hfg : f ≠ g
  htw : C.twin (f, s) = some (g, t)

namespace FlipAt

variable (D : C.FlipAt)

/-- Where an old side is after the flip: `j→k` moves to `(g, t)`, `i→l` to `(f, s)`. -/
def σ (x : F × Fin 3) : F × Fin 3 :=
  if x = (D.f, D.s + 1) then (D.g, D.t) else if x = (D.g, D.t + 1) then (D.f, D.s) else x

/-- The old side of a new side that is not the new edge. -/
def σi (y : F × Fin 3) : F × Fin 3 :=
  if y = (D.f, D.s) then (D.g, D.t + 1) else if y = (D.g, D.t) then (D.f, D.s + 1) else y

/-- Eq. 3, as `flipped_integers` evaluates it. -/
def nkl : ℤ :=
  flipNormal (C.n (D.f, D.s)) (C.n (D.f, D.s + 1)) (C.n (D.f, D.s + 2)) (C.n (D.g, D.t + 1))
    (C.n (D.g, D.t + 2))

/-- The side counts of the old face `ijk` from `i→j`. -/
def nA : Fin 3 → ℕ := fun k => C.np (D.f, D.s + k)
/-- The side counts of the old face `jil` from `j→i`. -/
def nB : Fin 3 → ℕ := fun k => C.np (D.g, D.t + k)
/-- `n⁻_ij`. -/
def m : ℕ := C.nm (D.f, D.s)

lemma ne_fg (a b : Fin 3) : (D.f, a) ≠ (D.g, b) := fun h => D.hfg (Prod.mk.inj h).1
lemma ne_gf (a b : Fin 3) : (D.g, a) ≠ (D.f, b) := fun h => D.hfg (Prod.mk.inj h).1.symm
lemma ne_ff {a b : Fin 3} (h : a ≠ b) : (D.f, a) ≠ (D.f, b) := fun e => h (Prod.mk.inj e).2
lemma ne_gg {a b : Fin 3} (h : a ≠ b) : (D.g, a) ≠ (D.g, b) := fun e => h (Prod.mk.inj e).2

lemma s1_ne (s : Fin 3) : s ≠ s + 1 := (fin3_add_one_ne s).symm
lemma s2_ne (s : Fin 3) : s ≠ s + 2 := (fin3_add_two_ne s).symm
lemma s12_ne (s : Fin 3) : s + 1 ≠ s + 2 := (fin3_add_two_ne_add_one s).symm

lemma σ_f1 : D.σ (D.f, D.s + 1) = (D.g, D.t) := by simp [σ]
lemma σ_g1 : D.σ (D.g, D.t + 1) = (D.f, D.s) := by simp [σ, D.ne_gf]
lemma σ_of {x : F × Fin 3} (h1 : x ≠ (D.f, D.s + 1)) (h2 : x ≠ (D.g, D.t + 1)) : D.σ x = x := by
  simp [σ, h1, h2]
lemma σi_f0 : D.σi (D.f, D.s) = (D.g, D.t + 1) := by simp [σi]
lemma σi_g0 : D.σi (D.g, D.t) = (D.f, D.s + 1) := by simp [σi, D.ne_gf]
lemma σi_of {y : F × Fin 3} (h1 : y ≠ (D.f, D.s)) (h2 : y ≠ (D.g, D.t)) : D.σi y = y := by
  simp [σi, h1, h2]

lemma σ_ne_f1 (x : F × Fin 3) : D.σ x ≠ (D.f, D.s + 1) := by
  unfold σ; split_ifs with h1 h2
  · exact D.ne_gf _ _
  · exact D.ne_ff (s1_ne _)
  · exact h1

lemma σ_ne_g1 (x : F × Fin 3) : D.σ x ≠ (D.g, D.t + 1) := by
  unfold σ; split_ifs with h1 h2
  · exact D.ne_gg (s1_ne _)
  · exact D.ne_fg _ _
  · exact h2

lemma σi_ne_f0 (y : F × Fin 3) : D.σi y ≠ (D.f, D.s) := by
  unfold σi; split_ifs with h1 h2
  · exact D.ne_gf _ _
  · exact D.ne_ff (s1_ne _).symm
  · exact h1

lemma σi_ne_g0 (y : F × Fin 3) : D.σi y ≠ (D.g, D.t) := by
  unfold σi; split_ifs with h1 h2
  · exact D.ne_gg (s1_ne _).symm
  · exact D.ne_fg _ _
  · exact h2

lemma σi_σ {x : F × Fin 3} (h1 : x ≠ (D.f, D.s)) (h2 : x ≠ (D.g, D.t)) : D.σi (D.σ x) = x := by
  by_cases e1 : x = (D.f, D.s + 1)
  · subst e1; rw [D.σ_f1, D.σi_g0]
  · by_cases e2 : x = (D.g, D.t + 1)
    · subst e2; rw [D.σ_g1, D.σi_f0]
    · rw [D.σ_of e1 e2, D.σi_of h1 h2]

lemma σ_σi {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    D.σ (D.σi y) = y := by
  by_cases e1 : y = (D.f, D.s)
  · subst e1; rw [D.σi_f0, D.σ_g1]
  · by_cases e2 : y = (D.g, D.t)
    · subst e2; rw [D.σi_g0, D.σ_f1]
    · rw [D.σi_of e1 e2, D.σ_of h1 h2]

/-- The twin of a side that is not on the flipped edge is not on it either. -/
lemma twin_ne_old (hC : C.Valid) {x z : F × Fin 3} (h : C.twin x = some z)
    (h1 : x ≠ (D.f, D.s)) (h2 : x ≠ (D.g, D.t)) : z ≠ (D.f, D.s) ∧ z ≠ (D.g, D.t) := by
  refine ⟨fun e => ?_, fun e => ?_⟩
  · subst e
    have := hC.twin_symm _ _ h
    rw [D.htw] at this
    exact h2 (Option.some.inj this).symm
  · subst e
    have h3 := hC.twin_symm _ _ h
    have h4 := hC.twin_symm _ _ D.htw
    rw [h3] at h4
    exact h1 (Option.some.inj h4)

end FlipAt

/-- **The flipped complex model** (`Complex::flip_edge`): new gluing and normal coordinates. -/
def flip (C : Cx F) (D : C.FlipAt) : Cx F where
  twin y := if y = (D.f, D.s + 1) then some (D.g, D.t + 1)
    else if y = (D.g, D.t + 1) then some (D.f, D.s + 1) else (C.twin (D.σi y)).map D.σ
  n y := if y = (D.f, D.s + 1) ∨ y = (D.g, D.t + 1) then D.nkl else C.n (D.σi y)

namespace FlipAt

variable (D : C.FlipAt)

lemma twin_f1 : (C.flip D).twin (D.f, D.s + 1) = some (D.g, D.t + 1) := by simp [flip]
lemma twin_g1 : (C.flip D).twin (D.g, D.t + 1) = some (D.f, D.s + 1) := by
  simp [flip, D.ne_gf]
lemma twin_of {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    (C.flip D).twin y = (C.twin (D.σi y)).map D.σ := by simp [flip, h1, h2]
lemma n_f1 : (C.flip D).n (D.f, D.s + 1) = D.nkl := by simp [flip]
lemma n_g1 : (C.flip D).n (D.g, D.t + 1) = D.nkl := by simp [flip]
lemma n_of {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    (C.flip D).n y = C.n (D.σi y) := by simp [flip, h1, h2]

end FlipAt

/-! ## Fin 3 bookkeeping -/

lemma fin3_rel (s k : Fin 3) : k = s ∨ k = s + 1 ∨ k = s + 2 := by
  fin_cases s <;> fin_cases k <;> decide

lemma fin3_sub_lits (s : Fin 3) : s + 1 - (s + 1) = 0 ∧ s + 2 - (s + 1) = 1 ∧ s - (s + 1) = 2 ∧
    s + 1 - s = 1 ∧ s + 2 - s = 2 ∧ s - s = 0 ∧ s + 0 = s := by
  fin_cases s <;> decide

/-- A valid complex model's faces are valid from any side on. -/
lemma Valid.face_at (hC : C.Valid) (f : F) (s : Fin 3) :
    faceValid (C.np (f, s)) (C.np (f, s + 1)) (C.np (f, s + 2)) = true := by
  have h := hC.face f
  rw [faceValid_nPlus] at h
  simp only [← np_cast] at h
  rcases fin3_cases s with rfl | rfl | rfl
  · exact h
  · rw [faceValid_rot] at h; exact h
  · rw [faceValid_rot, faceValid_rot] at h; exact h

namespace FlipAt

variable (D : C.FlipAt)

lemma nm_cast (hC : C.Valid) (x : F × Fin 3) : (C.n x : ℤ) = C.np x - C.nm x := by
  have := hC.n_ge x; simp only [np, nm]; split_ifs <;> omega

/-- The two old faces on the edge satisfy the hypotheses of the glued pair. -/
lemma hyp (hC : C.Valid) : Glue.Hyp D.nA D.nB D.m := by
  obtain ⟨-, -, -, -, -, -, l7⟩ := fin3_sub_lits D.s
  obtain ⟨-, -, -, -, -, -, l8⟩ := fin3_sub_lits D.t
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · simp only [nA, l7]; exact hC.face_at D.f D.s
  · simp only [nB, l8]; exact hC.face_at D.g D.t
  · simp only [nA, nB, l7, l8, np, hC.n_twin _ _ D.htw]
  · simp only [m, nm]; split_ifs <;> omega
  · intro h1; simp only [m, nm] at h1; simp only [nA, l7, np]; split_ifs at h1 with e
    · rw [e]; rfl

/-- **Eq. 3 is the count of the glued pair** (`Glue.Hyp.flipNormal_eq` on the complex). -/
lemma nkl_eq (hC : C.Valid) :
    D.nkl = (Glue.d D.nA D.nB D.m : ℤ) - Glue.Dkl D.nA D.nB := by
  obtain ⟨-, -, -, -, -, -, l7⟩ := fin3_sub_lits D.s
  obtain ⟨-, -, -, -, -, -, l8⟩ := fin3_sub_lits D.t
  apply (D.hyp hC).flipNormal_eq
  · simp only [nA, m, l7]; exact nm_cast hC _
  all_goals simp only [nA, nB, np_cast]

lemma Dkl_le (hC : C.Valid) (hn : -1 ≤ D.nkl) : Glue.Dkl D.nA D.nB ≤ 1 := by
  by_contra h
  have := (D.hyp hC).d_eq_zero (by omega)
  have e := D.nkl_eq hC
  rw [this] at e; omega

lemma np_nkl (hC : C.Valid) (hn : -1 ≤ D.nkl) : D.nkl.toNat = Glue.d D.nA D.nB D.m := by
  have e := D.nkl_eq hC
  have h1 := D.Dkl_le hC hn
  rcases Nat.lt_or_ge 0 (Glue.Dkl D.nA D.nB) with h | h
  · have := (D.hyp hC).d_eq_zero h; rw [this] at e ⊢; omega
  · omega

lemma nm_nkl (hC : C.Valid) (hn : -1 ≤ D.nkl) :
    (if D.nkl = -1 then 1 else 0) = Glue.Dkl D.nA D.nB := by
  have e := D.nkl_eq hC
  have h1 := D.Dkl_le hC hn
  rcases Nat.lt_or_ge 0 (Glue.Dkl D.nA D.nB) with h | h
  · have := (D.hyp hC).d_eq_zero h; rw [this] at e; split_ifs <;> omega
  · split_ifs <;> omega


lemma fn_A : C.fn D.f = fun k => D.nA (k - D.s) := by
  funext k; simp only [fn, nA]; congr 2; abel

lemma fn_B : C.fn D.g = fun k => D.nB (k - D.t) := by
  funext k; simp only [fn, nB]; congr 2; abel

/-- The side counts of the new face `f` = `F₁` rotated: slot `s + 1` is the new edge. -/
lemma fn_f (hC : C.Valid) (hn : -1 ≤ D.nkl) :
    (C.flip D).fn D.f = fun k => Glue.n1 D.nA D.nB D.m (k - (D.s + 1)) := by
  obtain ⟨l1, l2, l3, -, -, -, -⟩ := fin3_sub_lits D.s
  funext k
  rcases fin3_rel D.s k with rfl | rfl | rfl
  · simp only [fn, np, l3, Glue.n1, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    rw [D.n_of (D.ne_ff (s1_ne _)) (D.ne_fg _ _), D.σi_f0]; rfl
  · simp only [fn, np, l1, Glue.n1, Matrix.cons_val_zero, D.n_f1, D.np_nkl hC hn]
  · simp only [fn, np, l2, Glue.n1, Matrix.cons_val_one, Matrix.head_cons]
    rw [D.n_of (D.ne_ff (s12_ne _).symm) (D.ne_fg _ _),
      D.σi_of (D.ne_ff (s2_ne _).symm) (D.ne_fg _ _)]; rfl

/-- The side counts of the new face `g` = `F₂` rotated: slot `t + 1` is the new edge. -/
lemma fn_g (hC : C.Valid) (hn : -1 ≤ D.nkl) :
    (C.flip D).fn D.g = fun k => Glue.n2 D.nA D.nB D.m (k - (D.t + 1)) := by
  obtain ⟨l1, l2, l3, -, -, -, -⟩ := fin3_sub_lits D.t
  funext k
  rcases fin3_rel D.t k with rfl | rfl | rfl
  · simp only [fn, np, l3, Glue.n2, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    rw [D.n_of (D.ne_gf _ _) (D.ne_gg (s1_ne _)), D.σi_g0]; rfl
  · simp only [fn, np, l1, Glue.n2, Matrix.cons_val_zero, D.n_g1, D.np_nkl hC hn]
  · simp only [fn, np, l2, Glue.n2, Matrix.cons_val_one, Matrix.head_cons]
    rw [D.n_of (D.ne_gf _ _) (D.ne_gg (s12_ne _).symm),
      D.σi_of (D.ne_gf _ _) (D.ne_gg (s2_ne _).symm)]; rfl

lemma fn_other {h : F} (hf : h ≠ D.f) (hg : h ≠ D.g) : (C.flip D).fn h = C.fn h := by
  funext k
  have a : (h, k) ≠ (D.f, D.s + 1) := fun e => hf (Prod.mk.inj e).1
  have b : (h, k) ≠ (D.g, D.t + 1) := fun e => hg (Prod.mk.inj e).1
  have c : (h, k) ≠ (D.f, D.s) := fun e => hf (Prod.mk.inj e).1
  have d : (h, k) ≠ (D.g, D.t) := fun e => hg (Prod.mk.inj e).1
  simp only [fn, np, D.n_of a b, D.σi_of c d]

end FlipAt

lemma faceValid_shift (n : Fin 3 → ℕ) (r : Fin 3)
    (h : faceValid (n 0) (n 1) (n 2) = true) :
    faceValid (n (0 - r)) (n (1 - r)) (n (2 - r)) = true := by
  rcases fin3_cases r with rfl | rfl | rfl
  · simpa using h
  · rw [faceValid_rot, faceValid_rot] at h
    rw [show (0 : Fin 3) - 1 = 2 from rfl, show (1 : Fin 3) - 1 = 0 from rfl,
      show (2 : Fin 3) - 1 = 1 from rfl]
    exact h
  · rw [faceValid_rot] at h
    rw [show (0 : Fin 3) - 2 = 1 from rfl, show (1 : Fin 3) - 2 = 2 from rfl,
      show (2 : Fin 3) - 2 = 0 from rfl]
    exact h

lemma faceValid_of_fn {C' : Cx F} (f : F) (h : faceValid (C'.fn f 0) (C'.fn f 1) (C'.fn f 2) = true) :
    faceValid (C'.n (f, 0)) (C'.n (f, 1)) (C'.n (f, 2)) = true := by
  rw [faceValid_nPlus]; simpa only [fn, np_cast] using h

/-- **The flipped complex model is valid** (when `n_kl ≥ −1`, which the Rust code checks). The gluing
is an involution without fixed points, normal coordinates are `≥ −1`, equal on twins and `−1` on the
boundary, and the two new faces are valid: they are the faces `F₁`, `F₂` of the re-cut
(`Glue.Hyp.admissible1/2`). [formula: ic-flip-normal] -/
theorem flip_valid (hC : C.Valid) (D : C.FlipAt) (hn : -1 ≤ D.nkl) : (C.flip D).Valid := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · -- twin_symm
    intro y y' h
    by_cases e1 : y = (D.f, D.s + 1)
    · subst e1; rw [D.twin_f1] at h; cases h; exact D.twin_g1
    by_cases e2 : y = (D.g, D.t + 1)
    · subst e2; rw [D.twin_g1] at h; cases h; exact D.twin_f1
    rw [D.twin_of e1 e2] at h
    obtain ⟨z, hz, rfl⟩ := Option.map_eq_some_iff.mp h
    obtain ⟨z1, z2⟩ := D.twin_ne_old hC hz (D.σi_ne_f0 y) (D.σi_ne_g0 y)
    rw [D.twin_of (D.σ_ne_f1 z) (D.σ_ne_g1 z), D.σi_σ z1 z2, hC.twin_symm _ _ hz]
    simp [D.σ_σi e1 e2]
  · -- twin_ne
    intro y h
    by_cases e1 : y = (D.f, D.s + 1)
    · subst e1; rw [D.twin_f1] at h; exact D.ne_gf _ _ (Option.some.inj h)
    by_cases e2 : y = (D.g, D.t + 1)
    · subst e2; rw [D.twin_g1] at h; exact D.ne_fg _ _ (Option.some.inj h)
    rw [D.twin_of e1 e2] at h
    obtain ⟨z, hz, hzy⟩ := Option.map_eq_some_iff.mp h
    obtain ⟨z1, z2⟩ := D.twin_ne_old hC hz (D.σi_ne_f0 y) (D.σi_ne_g0 y)
    have : z = D.σi y := by rw [← hzy, D.σi_σ z1 z2]
    rw [this] at hz
    exact hC.twin_ne _ hz
  · -- n_ge
    intro y
    by_cases e : y = (D.f, D.s + 1) ∨ y = (D.g, D.t + 1)
    · simp only [flip, e, ite_true]; exact hn
    · simp only [flip, e, ite_false]; exact hC.n_ge _
  · -- n_twin
    intro y y' h
    by_cases e1 : y = (D.f, D.s + 1)
    · subst e1; rw [D.twin_f1] at h; cases h; rw [D.n_f1, D.n_g1]
    by_cases e2 : y = (D.g, D.t + 1)
    · subst e2; rw [D.twin_g1] at h; cases h; rw [D.n_f1, D.n_g1]
    rw [D.twin_of e1 e2] at h
    obtain ⟨z, hz, rfl⟩ := Option.map_eq_some_iff.mp h
    obtain ⟨z1, z2⟩ := D.twin_ne_old hC hz (D.σi_ne_f0 y) (D.σi_ne_g0 y)
    rw [D.n_of (D.σ_ne_f1 z) (D.σ_ne_g1 z), D.σi_σ z1 z2, D.n_of e1 e2, hC.n_twin _ _ hz]
  · -- n_bdry
    intro y h
    by_cases e1 : y = (D.f, D.s + 1)
    · subst e1; rw [D.twin_f1] at h; cases h
    by_cases e2 : y = (D.g, D.t + 1)
    · subst e2; rw [D.twin_g1] at h; cases h
    rw [D.twin_of e1 e2, Option.map_eq_none_iff] at h
    rw [D.n_of e1 e2]; exact hC.n_bdry _ h
  · -- the faces
    intro h
    have hD := D.Dkl_le hC hn
    by_cases hf : h = D.f
    · subst hf
      apply faceValid_of_fn
      rw [D.fn_f hC hn]
      have := (FaceModel.exists_valid_iff _).1 ⟨_, rfl, canon_valid ((D.hyp hC).admissible1 hD)⟩
      exact faceValid_shift _ _ this
    by_cases hg : h = D.g
    · subst hg
      apply faceValid_of_fn
      rw [D.fn_g hC hn]
      have := (FaceModel.exists_valid_iff _).1 ⟨_, rfl, canon_valid ((D.hyp hC).admissible2 hD)⟩
      exact faceValid_shift _ _ this
    apply faceValid_of_fn
    rw [D.fn_other hf hg]
    have := hC.face h
    rw [faceValid_nPlus] at this
    simpa only [fn, np_cast] using this


/-! ## The corners after the flip -/

namespace FlipAt

variable (D : C.FlipAt)

/-- The old corner at the vertex of a new one (the halfedge's tail): the new edge `l→k` starts at
`l` (old corner `(g, t + 2)`), `k→l` at `k` (old corner `(f, s + 2)`); every other side is the
old side `σi y`. -/
def oc (y : F × Fin 3) : F × Fin 3 :=
  if y = (D.f, D.s + 1) then (D.g, D.t + 2) else if y = (D.g, D.t + 1) then (D.f, D.s + 2)
  else D.σi y

/-- Eq. 4 at `k`, as `updated_roundabout` computes it (from `r_ki`, the emanating count of corner
`k` of `kil` from the new counts, `n⁻_ki`, the modulus). -/
def rk (L : Labels F) : ℕ :=
  (roundaboutUpdate (L.r (D.f, D.s + 2)) (emanating (C.n (D.g, D.t + 1)) (C.n (D.f, D.s + 2)) D.nkl)
    (nMinus (C.n (D.f, D.s + 2))) (L.md (D.f, D.s + 2))).toNat

/-- Eq. 4 at `l` (from `r_lj`, corner `l` of `ljk`, `n⁻_lj`). -/
def rl (L : Labels F) : ℕ :=
  (roundaboutUpdate (L.r (D.g, D.t + 2)) (emanating (C.n (D.f, D.s + 1)) (C.n (D.g, D.t + 2)) D.nkl)
    (nMinus (C.n (D.g, D.t + 2))) (L.md (D.g, D.t + 2))).toNat

/-- The labels after the flip: the same input vertices, degrees and moduli; `r_kl`, `r_lk` by
Eq. 4; every other halfedge keeps its roundabout (the two that change faces move with it). -/
def labels (L : Labels F) : Labels F where
  input y := L.input (D.oc y)
  deg y := L.deg (D.oc y)
  md y := L.md (D.oc y)
  r y := if y = (D.g, D.t + 1) then D.rk L else if y = (D.f, D.s + 1) then D.rl L else L.r (D.σi y)

lemma oc_f1 : D.oc (D.f, D.s + 1) = (D.g, D.t + 2) := by simp [oc]
lemma oc_g1 : D.oc (D.g, D.t + 1) = (D.f, D.s + 2) := by simp [oc, D.ne_gf]
lemma oc_of {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    D.oc y = D.σi y := by simp [oc, h1, h2]

lemma r_g1 (L : Labels F) : (D.labels L).r (D.g, D.t + 1) = D.rk L := by simp [labels]
lemma r_f1 (L : Labels F) : (D.labels L).r (D.f, D.s + 1) = D.rl L := by simp [labels, D.ne_fg]
lemma r_of (L : Labels F) {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    (D.labels L).r y = L.r (D.σi y) := by simp [labels, h1, h2]

/-- The labels at a moved or kept side are the old ones. -/
lemma labels_σ (L : Labels F) {x : F × Fin 3} (h1 : x ≠ (D.f, D.s)) (h2 : x ≠ (D.g, D.t)) :
    ((D.labels L).input (D.σ x) ↔ L.input x) ∧ (D.labels L).deg (D.σ x) = L.deg x ∧
      (D.labels L).md (D.σ x) = L.md x ∧ (D.labels L).r (D.σ x) = L.r x := by
  have a := D.σ_ne_f1 x
  have b := D.σ_ne_g1 x
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp only [labels, D.oc_of a b, D.σi_σ h1 h2]
  · simp only [labels, D.oc_of a b, D.σi_σ h1 h2]
  · simp only [labels, D.oc_of a b, D.σi_σ h1 h2]
  · rw [D.r_of L a b, D.σi_σ h1 h2]

/-! ### Rotation around the vertices after the flip -/

lemma rot_old_g1 (hC : C.Valid) : C.rot (D.g, D.t + 1) = some (D.f, D.s) := by
  simp only [rot, fin3_add_one_add_two]; exact hC.twin_symm _ _ D.htw

lemma rot_old_f1 : C.rot (D.f, D.s + 1) = some (D.g, D.t) := by
  simp only [rot, fin3_add_one_add_two]; exact D.htw

lemma rot'_f0 : (C.flip D).rot (D.f, D.s) = (C.rot (D.f, D.s)).map D.σ := by
  simp only [rot]
  rw [D.twin_of (D.ne_ff (s12_ne _).symm) (D.ne_fg _ _),
    D.σi_of (D.ne_ff (s2_ne _).symm) (D.ne_fg _ _)]

lemma rot'_f1 : (C.flip D).rot (D.f, D.s + 1) = (C.rot (D.g, D.t + 2)).map D.σ := by
  simp only [rot, fin3_add_one_add_two, fin3_add_two_add_two]
  rw [D.twin_of (D.ne_ff (s1_ne _)) (D.ne_fg _ _), D.σi_f0]

lemma rot'_f2 : (C.flip D).rot (D.f, D.s + 2) = some (D.g, D.t + 1) := by
  simp only [rot, fin3_add_two_add_two]; exact D.twin_f1

lemma rot'_g0 : (C.flip D).rot (D.g, D.t) = (C.rot (D.g, D.t)).map D.σ := by
  simp only [rot]
  rw [D.twin_of (D.ne_gf _ _) (D.ne_gg (s12_ne _).symm),
    D.σi_of (D.ne_gf _ _) (D.ne_gg (s2_ne _).symm)]

lemma rot'_g1 : (C.flip D).rot (D.g, D.t + 1) = (C.rot (D.f, D.s + 2)).map D.σ := by
  simp only [rot, fin3_add_one_add_two, fin3_add_two_add_two]
  rw [D.twin_of (D.ne_gf _ _) (D.ne_gg (s1_ne _)), D.σi_g0]

lemma rot'_g2 : (C.flip D).rot (D.g, D.t + 2) = some (D.f, D.s + 1) := by
  simp only [rot, fin3_add_two_add_two]; exact D.twin_g1

lemma rot'_O {c : F × Fin 3} (hf : c.1 ≠ D.f) (hg : c.1 ≠ D.g) :
    (C.flip D).rot c = (C.rot c).map D.σ := by
  simp only [rot]
  rw [D.twin_of (fun e => hf (Prod.mk.inj e).1) (fun e => hg (Prod.mk.inj e).1),
    D.σi_of (fun e => hf (Prod.mk.inj e).1) (fun e => hg (Prod.mk.inj e).1)]

lemma twin'_f0 : (C.flip D).twin (D.f, D.s) = (C.twin (D.g, D.t + 1)).map D.σ := by
  rw [D.twin_of (D.ne_ff (s1_ne _)) (D.ne_fg _ _), D.σi_f0]
lemma twin'_f2 : (C.flip D).twin (D.f, D.s + 2) = (C.twin (D.f, D.s + 2)).map D.σ := by
  rw [D.twin_of (D.ne_ff (s12_ne _).symm) (D.ne_fg _ _),
    D.σi_of (D.ne_ff (s2_ne _).symm) (D.ne_fg _ _)]
lemma twin'_g0 : (C.flip D).twin (D.g, D.t) = (C.twin (D.f, D.s + 1)).map D.σ := by
  rw [D.twin_of (D.ne_gf _ _) (D.ne_gg (s1_ne _)), D.σi_g0]
lemma twin'_g2 : (C.flip D).twin (D.g, D.t + 2) = (C.twin (D.g, D.t + 2)).map D.σ := by
  rw [D.twin_of (D.ne_gf _ _) (D.ne_gg (s12_ne _).symm),
    D.σi_of (D.ne_gf _ _) (D.ne_gg (s2_ne _).symm)]
lemma twin'_O {c : F × Fin 3} (hf : c.1 ≠ D.f) (hg : c.1 ≠ D.g) :
    (C.flip D).twin c = (C.twin c).map D.σ := by
  rw [D.twin_of (fun e => hf (Prod.mk.inj e).1) (fun e => hg (Prod.mk.inj e).1),
    D.σi_of (fun e => hf (Prod.mk.inj e).1) (fun e => hg (Prod.mk.inj e).1)]

/-- The old rotation of these corners never lands on the flipped edge. -/
lemma rot_ne_old (hC : C.Valid) {c x : F × Fin 3} (h : C.rot c = some x)
    (h1 : (c.1, c.2 + 2) ≠ (D.f, D.s)) (h2 : (c.1, c.2 + 2) ≠ (D.g, D.t)) :
    x ≠ (D.f, D.s) ∧ x ≠ (D.g, D.t) := D.twin_ne_old hC h h1 h2

end FlipAt


/-! ### Counts at the corners after the flip -/

lemma roundaboutUpdate_nat (r e n md : ℕ) :
    (roundaboutUpdate r e n md).toNat = (r + e + n) % md := by
  unfold roundaboutUpdate
  have e : (r : ℤ) + e + n = ((r + e + n : ℕ) : ℤ) := by push_cast; ring
  split_ifs with h
  · simp only [beq_iff_eq] at h; subst h; rw [e, Int.toNat_natCast, Nat.mod_zero]
  · rw [e, ← Int.natCast_mod, Int.toNat_natCast]

namespace FlipAt

variable (D : C.FlipAt)

lemma emc_A (hC : C.Valid) (k : Fin 3) : C.emc (D.f, D.s + k) = Glue.ea D.nA k := by
  rw [hC.emc_eq]
  simp only
  rw [D.fn_A, show (fun k => D.nA (k - D.s)) = fun k => D.nA (k + -D.s) by
    funext k; rw [sub_eq_add_neg], em_rot]
  simp only [Glue.ea]; congr 1; abel

lemma emc_B (hC : C.Valid) (k : Fin 3) : C.emc (D.g, D.t + k) = Glue.eb D.nB k := by
  rw [hC.emc_eq]
  simp only
  rw [D.fn_B, show (fun k => D.nB (k - D.t)) = fun k => D.nB (k + -D.t) by
    funext k; rw [sub_eq_add_neg], em_rot]
  simp only [Glue.eb]; congr 1; abel

lemma emc'_f (hC : C.Valid) (hn : -1 ≤ D.nkl) (k : Fin 3) :
    (C.flip D).emc (D.f, k) = em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) (k - (D.s + 1)) := by
  rw [(flip_valid hC D hn).emc_eq]
  simp only
  rw [D.fn_f hC hn, show (fun k => Glue.n1 D.nA D.nB D.m (k - (D.s + 1))) =
    fun k => Glue.n1 D.nA D.nB D.m (k + -(D.s + 1)) by funext k; rw [sub_eq_add_neg], em_rot,
    (D.hyp hC).cc1 (D.Dkl_le hC hn), sub_eq_add_neg]

lemma emc'_g (hC : C.Valid) (hn : -1 ≤ D.nkl) (k : Fin 3) :
    (C.flip D).emc (D.g, k) = em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) (k - (D.t + 1)) := by
  rw [(flip_valid hC D hn).emc_eq]
  simp only
  rw [D.fn_g hC hn, show (fun k => Glue.n2 D.nA D.nB D.m (k - (D.t + 1))) =
    fun k => Glue.n2 D.nA D.nB D.m (k + -(D.t + 1)) by funext k; rw [sub_eq_add_neg], em_rot,
    (D.hyp hC).cc2 (D.Dkl_le hC hn), sub_eq_add_neg]

lemma emc'_O (hC : C.Valid) (hn : -1 ≤ D.nkl) {c : F × Fin 3} (hf : c.1 ≠ D.f) (hg : c.1 ≠ D.g) :
    (C.flip D).emc c = C.emc c := by
  rw [(flip_valid hC D hn).emc_eq, hC.emc_eq, D.fn_other hf hg]

lemma nm'_new (hC : C.Valid) (hn : -1 ≤ D.nkl) :
    (C.flip D).nm (D.f, D.s + 1) = Glue.Dkl D.nA D.nB ∧
      (C.flip D).nm (D.g, D.t + 1) = Glue.Dkl D.nA D.nB := by
  have := D.nm_nkl hC hn
  simp only [nm, D.n_f1, D.n_g1]; exact ⟨this, this⟩

lemma nm'_of {y : F × Fin 3} (h1 : y ≠ (D.f, D.s + 1)) (h2 : y ≠ (D.g, D.t + 1)) :
    (C.flip D).nm y = C.nm (D.σi y) := by
  simp only [nm, D.n_of h1 h2]

lemma rk_eq (hC : C.Valid) (hn : -1 ≤ D.nkl) (L : Labels F) :
    D.rk L = (L.r (D.f, D.s + 2) + em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) 1 +
      C.nm (D.f, D.s + 2)) % L.md (D.f, D.s + 2) := by
  have e := (flip_valid hC D hn).emc_eq_emanating (D.f, D.s + 2)
  simp only [fin3_add_two_add_one, fin3_add_two_add_two] at e
  rw [D.n_of (D.ne_ff (s1_ne _)) (D.ne_fg _ _), D.σi_f0,
    D.n_of (D.ne_ff (s12_ne _).symm) (D.ne_fg _ _), D.σi_of (D.ne_ff (s2_ne _).symm) (D.ne_fg _ _),
    D.n_f1, D.emc'_f hC hn] at e
  rw [show D.s + 2 - (D.s + 1) = 1 from (fin3_sub_lits D.s).2.1] at e
  unfold rk
  rw [← e, ← nm_eq hC, roundaboutUpdate_nat]

lemma rl_eq (hC : C.Valid) (hn : -1 ≤ D.nkl) (L : Labels F) :
    D.rl L = (L.r (D.g, D.t + 2) + em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) 1 +
      C.nm (D.g, D.t + 2)) % L.md (D.g, D.t + 2) := by
  have e := (flip_valid hC D hn).emc_eq_emanating (D.g, D.t + 2)
  simp only [fin3_add_two_add_one, fin3_add_two_add_two] at e
  rw [D.n_of (D.ne_gf _ _) (D.ne_gg (s1_ne _)), D.σi_g0,
    D.n_of (D.ne_gf _ _) (D.ne_gg (s12_ne _).symm), D.σi_of (D.ne_gf _ _) (D.ne_gg (s2_ne _).symm),
    D.n_g1, D.emc'_g hC hn] at e
  rw [show D.t + 2 - (D.t + 1) = 1 from (fin3_sub_lits D.t).2.1] at e
  unfold rl
  rw [← e, ← nm_eq hC, roundaboutUpdate_nat]

end FlipAt


namespace FlipAt

variable (D : C.FlipAt)

lemma corner_cases (c : F × Fin 3) :
    c = (D.f, D.s) ∨ c = (D.f, D.s + 1) ∨ c = (D.f, D.s + 2) ∨ c = (D.g, D.t) ∨
      c = (D.g, D.t + 1) ∨ c = (D.g, D.t + 2) ∨ (c.1 ≠ D.f ∧ c.1 ≠ D.g) := by
  obtain ⟨h, k⟩ := c
  by_cases hf : h = D.f
  · subst hf; rcases fin3_rel D.s k with rfl | rfl | rfl <;> simp
  by_cases hg : h = D.g
  · subst hg; rcases fin3_rel D.t k with rfl | rfl | rfl <;> simp
  exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr ⟨hf, hg⟩)))))

end FlipAt

/-- **`ic-roundabout-update`: Eq. 4 keeps the roundabouts consistent.** If the roundabouts of a valid
state agree with the counter-clockwise order of the curve ends (`RoundaboutsOK`, what `certify`
checks; `RoundaboutsOK.lab_next`: the ends get consecutive local indices), then after the flip of
an edge with `n_kl ≥ −1` they do too, with `r_kl` and `r_lk` from Eq. 4 (`FlipAt.rk`, `FlipAt.rl`:
the oracle's `roundaboutUpdate` of `r_ki` (`r_lj`), the emanating count of corner `k` of `kil`
(`l` of `ljk`) from the new normal coordinates, and `n⁻_ki` (`n⁻_lj`), exactly the arguments of
the Rust `updated_roundabout`) and every other halfedge keeping its roundabout (the halfedges
`i→l`, `j→k` moving with their data). Since the curves are unchanged (`Glue.Hyp.outer`) and the
labels of the ends around each vertex stay consecutive, Eq. 4 is the local index of the first input
halfedge at or after `k→l`. [formula: ic-roundabout-update] -/
theorem flip_roundabouts (hC : C.Valid) (D : C.FlipAt) (hn : -1 ≤ D.nkl) {L : Labels F}
    (hL : C.RoundaboutsOK L) : (C.flip D).RoundaboutsOK (D.labels L) := by
  have hD := D.Dkl_le hC hn
  have H := D.hyp hC
  obtain ⟨l1, l2, -, -, -, -, l7⟩ := fin3_sub_lits D.s
  obtain ⟨l1', l2', -, -, -, -, l7'⟩ := fin3_sub_lits D.t
  -- the old corners and their counts
  have R1 := D.rot_old_g1 hC
  have R2 := D.rot_old_f1
  have eA0 := D.emc_A hC 0
  have eA1 := D.emc_A hC 1
  have eA2 := D.emc_A hC 2
  have eB0 := D.emc_B hC 0
  have eB1 := D.emc_B hC 1
  have eB2 := D.emc_B hC 2
  rw [l7] at eA0
  rw [l7'] at eB0
  have Ek := H.em_k hD
  have El := H.em_l hD
  have Ei := H.em_i hD
  have Ej := H.em_j hD
  -- the new corners' counts
  have f0 := D.emc'_f hC hn D.s
  have f1 := D.emc'_f hC hn (D.s + 1)
  have f2 := D.emc'_f hC hn (D.s + 2)
  have g0 := D.emc'_g hC hn D.t
  have g1 := D.emc'_g hC hn (D.t + 1)
  have g2 := D.emc'_g hC hn (D.t + 2)
  obtain ⟨-, -, s3, -, -, -, -⟩ := fin3_sub_lits D.s
  obtain ⟨-, -, t3, -, -, -, -⟩ := fin3_sub_lits D.t
  rw [s3] at f0; rw [l1] at f1; rw [l2] at f2
  rw [t3] at g0; rw [l1'] at g1; rw [l2'] at g2
  obtain ⟨nf1, ng1⟩ := D.nm'_new hC hn
  have nf0 := D.nm'_of (D.ne_ff (FlipAt.s1_ne D.s)) (D.ne_fg D.s (D.t + 1))
  have nf2 := D.nm'_of (D.ne_ff (FlipAt.s12_ne D.s).symm) (D.ne_fg (D.s + 2) (D.t + 1))
  have ng0 := D.nm'_of (D.ne_gf D.t (D.s + 1)) (D.ne_gg (FlipAt.s1_ne D.t))
  have ng2 := D.nm'_of (D.ne_gf (D.t + 2) (D.s + 1)) (D.ne_gg (FlipAt.s12_ne D.t).symm)
  rw [D.σi_f0] at nf0
  rw [D.σi_of (D.ne_ff (FlipAt.s2_ne D.s).symm) (D.ne_fg _ _)] at nf2
  rw [D.σi_g0] at ng0
  rw [D.σi_of (D.ne_gf _ _) (D.ne_gg (FlipAt.s2_ne D.t).symm)] at ng2
  have nm_gt : C.nm (D.g, D.t) = D.m := by
    simp only [nm, FlipAt.m, hC.n_twin _ _ D.htw]
  have hrk := D.rk_eq hC hn L
  have hrl := D.rl_eq hC hn L
  -- the weights `n⁻ + e` of the corners, after and before
  have W0 : (C.flip D).w (D.f, D.s) = C.w (D.g, D.t + 1) + C.w (D.f, D.s) := by
    simp only [w]
    rw [nf0, f0, Ei, eA0, eB1]
    simp only [FlipAt.m]
    ring
  have W1 : (C.flip D).w (D.f, D.s + 1) =
      Glue.Dkl D.nA D.nB + em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) 0 := by
    simp only [w, nf1, f1]
  have W2 : (C.flip D).w (D.f, D.s + 2) =
      C.nm (D.f, D.s + 2) + em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) 1 := by
    simp only [w, nf2, f2]
  have V0 : (C.flip D).w (D.g, D.t) = C.w (D.f, D.s + 1) + C.w (D.g, D.t) := by
    simp only [w]
    rw [ng0, g0, Ej, eA1, eB0, nm_gt]
    ring
  have V1 : (C.flip D).w (D.g, D.t + 1) =
      Glue.Dkl D.nA D.nB + em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) 0 := by
    simp only [w, ng1, g1]
  have V2 : (C.flip D).w (D.g, D.t + 2) =
      C.nm (D.g, D.t + 2) + em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) 1 := by
    simp only [w, ng2, g2]
  have Wk : C.w (D.f, D.s + 2) = C.nm (D.f, D.s + 2) +
      em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) 1 + Glue.Dkl D.nA D.nB +
      em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) 0 := by
    simp only [w]; rw [eA2, Ek]; ring
  have Wl : C.w (D.g, D.t + 2) = C.nm (D.g, D.t + 2) +
      em (Glue.n1 D.nA D.nB D.m) (Glue.c1 D.nA D.nB) 0 + Glue.Dkl D.nA D.nB +
      em (Glue.n2 D.nA D.nB D.m) (Glue.c2 D.nA D.nB) 1 := by
    simp only [w]; rw [eB2, El]; ring
  have nq : ∀ c : F × Fin 3, c.1 ≠ D.f → c.1 ≠ D.g →
      c ≠ (D.f, D.s + 1) ∧ c ≠ (D.g, D.t + 1) ∧ c ≠ (D.f, D.s) ∧ c ≠ (D.g, D.t) := by
    intro c hf hg
    refine ⟨fun e => hf ?_, fun e => hg ?_, fun e => hf ?_, fun e => hg ?_⟩ <;> rw [e]
  have WO : ∀ c : F × Fin 3, c.1 ≠ D.f → c.1 ≠ D.g → (C.flip D).w c = C.w c := by
    intro c hf hg
    obtain ⟨a, b, a', b'⟩ := nq c hf hg
    simp only [w, D.nm'_of a b, D.emc'_O hC hn hf hg, D.σi_of a' b']
  -- the labels at the corners
  have oc0 : D.oc (D.f, D.s) = (D.g, D.t + 1) := by
    rw [D.oc_of (D.ne_ff (FlipAt.s1_ne _)) (D.ne_fg _ _), D.σi_f0]
  have oc2 : D.oc (D.f, D.s + 2) = (D.f, D.s + 2) := by
    rw [D.oc_of (D.ne_ff (FlipAt.s12_ne _).symm) (D.ne_fg _ _),
      D.σi_of (D.ne_ff (FlipAt.s2_ne _).symm) (D.ne_fg _ _)]
  have od0 : D.oc (D.g, D.t) = (D.f, D.s + 1) := by
    rw [D.oc_of (D.ne_gf _ _) (D.ne_gg (FlipAt.s1_ne _)), D.σi_g0]
  have od2 : D.oc (D.g, D.t + 2) = (D.g, D.t + 2) := by
    rw [D.oc_of (D.ne_gf _ _) (D.ne_gg (FlipAt.s12_ne _).symm),
      D.σi_of (D.ne_gf _ _) (D.ne_gg (FlipAt.s2_ne _).symm)]
  have oO : ∀ c : F × Fin 3, c.1 ≠ D.f → c.1 ≠ D.g → D.oc c = c := by
    intro c hf hg
    obtain ⟨a, b, a', b'⟩ := nq c hf hg
    rw [D.oc_of a b, D.σi_of a' b']
  have r0 : (D.labels L).r (D.f, D.s) = L.r (D.g, D.t + 1) := by
    rw [D.r_of L (D.ne_ff (FlipAt.s1_ne _)) (D.ne_fg _ _), D.σi_f0]
  have r2 : (D.labels L).r (D.f, D.s + 2) = L.r (D.f, D.s + 2) := by
    rw [D.r_of L (D.ne_ff (FlipAt.s12_ne _).symm) (D.ne_fg _ _),
      D.σi_of (D.ne_ff (FlipAt.s2_ne _).symm) (D.ne_fg _ _)]
  have q0 : (D.labels L).r (D.g, D.t) = L.r (D.f, D.s + 1) := by
    rw [D.r_of L (D.ne_gf _ _) (D.ne_gg (FlipAt.s1_ne _)), D.σi_g0]
  have q2 : (D.labels L).r (D.g, D.t + 2) = L.r (D.g, D.t + 2) := by
    rw [D.r_of L (D.ne_gf _ _) (D.ne_gg (FlipAt.s12_ne _).symm),
      D.σi_of (D.ne_gf _ _) (D.ne_gg (FlipAt.s2_ne _).symm)]
  have rO : ∀ c : F × Fin 3, c.1 ≠ D.f → c.1 ≠ D.g → (D.labels L).r c = L.r c := by
    intro c hf hg
    obtain ⟨a, b, a', b'⟩ := nq c hf hg
    rw [D.r_of L a b, D.σi_of a' b']
  have rf1 := D.r_f1 L
  have rg1 := D.r_g1 L
  -- the rotation after the flip, and where the old one lands
  have nold : ∀ {c x : F × Fin 3}, C.rot c = some x → (c.1, c.2 + 2) ≠ (D.f, D.s) →
      (c.1, c.2 + 2) ≠ (D.g, D.t) → x ≠ (D.f, D.s) ∧ x ≠ (D.g, D.t) :=
    fun h a b => D.rot_ne_old hC h a b
  have lσ := fun {x : F × Fin 3} (h1 : x ≠ (D.f, D.s)) (h2 : x ≠ (D.g, D.t)) => D.labels_σ L h1 h2
  -- vertex invariance of the labels
  have vtx : ∀ {α : Type} (X : F × Fin 3 → α), (∀ u v, C.rot u = some v → X v = X u) →
      ∀ c c', (C.flip D).rot c = some c' → X (D.oc c') = X (D.oc c) := by
    intro α X hX c c' h
    have ocσ : ∀ {x : F × Fin 3}, x ≠ (D.f, D.s) → x ≠ (D.g, D.t) → D.oc (D.σ x) = x :=
      fun h1 h2 => by rw [D.oc_of (D.σ_ne_f1 _) (D.σ_ne_g1 _), D.σi_σ h1 h2]
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [D.rot'_f0] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_ff (FlipAt.s2_ne _).symm) (D.ne_fg _ _)
      rw [ocσ x1 x2, oc0, hX _ _ hx, hX _ _ R1]
    · rw [D.rot'_f1] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_gf _ _) (by simp only [fin3_add_two_add_two]; exact D.ne_gg (FlipAt.s1_ne _).symm)
      rw [ocσ x1 x2, D.oc_f1, hX _ _ hx]
    · rw [D.rot'_f2] at h; cases h; rw [D.oc_g1, oc2]
    · rw [D.rot'_g0] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_gf _ _) (D.ne_gg (FlipAt.s2_ne _).symm)
      rw [ocσ x1 x2, od0, hX _ _ hx, hX _ _ R2]
    · rw [D.rot'_g1] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (by simp only [fin3_add_two_add_two]; exact D.ne_ff (FlipAt.s1_ne _).symm) (D.ne_fg _ _)
      rw [ocσ x1 x2, D.oc_g1, hX _ _ hx]
    · rw [D.rot'_g2] at h; cases h; rw [D.oc_f1, od2]
    · rw [D.rot'_O hf hg] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (fun e => hf (by rw [(Prod.mk.inj e).1])) (fun e => hg (by rw [(Prod.mk.inj e).1]))
      rw [ocσ x1 x2, oO c hf hg, hX _ _ hx]
  have inp0 : ∀ hc : L.input (D.g, D.t + 1), L.input (D.f, D.s) := fun hc =>
    (hL.input_rot _ _ R1).2 hc
  have inp1 : ∀ hc : L.input (D.f, D.s + 1), L.input (D.g, D.t) := fun hc =>
    (hL.input_rot _ _ R2).2 hc
  have md0 := hL.md_rot _ _ R1
  have md1 := hL.md_rot _ _ R2
  have dg0 := hL.deg_rot _ _ R1
  have dg1 := hL.deg_rot _ _ R2
  -- a bound used for Eq. 4 at a boundary vertex
  have bound : ∀ x : F × Fin 3, L.input x → L.md x = 0 → L.r x + C.w x < L.deg x := by
    intro x hx hm
    rcases e : C.rot x with _ | y
    · have := hL.r_last x hx e; omega
    · have h1 := hL.r_next x y hx e
      have h2 := hL.r_lt y ((hL.input_rot x y e).2 hx)
      rw [hm, Nat.mod_zero] at h1
      rw [hL.deg_rot x y e] at h2; omega
  refine
    { input_rot := fun c c' h => ?_
      deg_rot := fun c c' h => vtx (fun y => L.deg y) (fun u v huv => hL.deg_rot u v huv) c c' h
      md_rot := fun c c' h => vtx (fun y => L.md y) (fun u v huv => hL.md_rot u v huv) c c' h
      md_eq := fun c hc => hL.md_eq _ hc
      md_bdry := fun c hc h => ?_
      r_lt := fun c hc => ?_
      r_next := fun c c' hc h => ?_
      r_last := fun c hc h => ?_
      r_first := fun c hc h => ?_ }
  · exact Iff.of_eq (vtx (fun y => L.input y) (fun u v huv => propext (hL.input_rot u v huv)) c c' h)
  · -- md_bdry
    show L.md (D.oc c) = 0
    change L.input (D.oc c) at hc
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [D.rot'_f0, Option.map_eq_none_iff] at h
      rw [oc0] at hc ⊢; rw [← md0]; exact hL.md_bdry _ (inp0 hc) h
    · rw [D.rot'_f1, Option.map_eq_none_iff] at h
      rw [D.oc_f1] at hc ⊢; exact hL.md_bdry _ hc h
    · rw [D.rot'_f2] at h; cases h
    · rw [D.rot'_g0, Option.map_eq_none_iff] at h
      rw [od0] at hc ⊢; rw [← md1]; exact hL.md_bdry _ (inp1 hc) h
    · rw [D.rot'_g1, Option.map_eq_none_iff] at h
      rw [D.oc_g1] at hc ⊢; exact hL.md_bdry _ hc h
    · rw [D.rot'_g2] at h; cases h
    · rw [D.rot'_O hf hg, Option.map_eq_none_iff] at h
      rw [oO c hf hg] at hc ⊢; exact hL.md_bdry _ hc h
  · -- r_lt
    show (D.labels L).r c < L.deg (D.oc c)
    change L.input (D.oc c) at hc
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [r0, oc0]; rw [oc0] at hc; exact hL.r_lt _ hc
    · rw [rf1, D.oc_f1, hrl]; rw [D.oc_f1] at hc
      have hr := hL.r_lt _ hc
      rcases hL.md_eq _ hc with h0 | h0
      · rw [h0, Nat.mod_zero]; have := bound _ hc h0; rw [Wl] at this; omega
      · rw [h0]; exact Nat.mod_lt _ (by omega)
    · rw [r2, oc2]; rw [oc2] at hc; exact hL.r_lt _ hc
    · rw [q0, od0]; rw [od0] at hc; exact hL.r_lt _ hc
    · rw [rg1, D.oc_g1, hrk]; rw [D.oc_g1] at hc
      have hr := hL.r_lt _ hc
      rcases hL.md_eq _ hc with h0 | h0
      · rw [h0, Nat.mod_zero]; have := bound _ hc h0; rw [Wk] at this; omega
      · rw [h0]; exact Nat.mod_lt _ (by omega)
    · rw [q2, od2]; rw [od2] at hc; exact hL.r_lt _ hc
    · rw [rO c hf hg, oO c hf hg]; rw [oO c hf hg] at hc; exact hL.r_lt _ hc
  · -- r_next
    show ((D.labels L).r c + (C.flip D).w c) % L.md (D.oc c) = (D.labels L).r c'
    change L.input (D.oc c) at hc
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [D.rot'_f0] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_ff (FlipAt.s2_ne _).symm) (D.ne_fg _ _)
      rw [(lσ x1 x2).2.2.2, r0, oc0, W0, ← hL.r_next _ _ (inp0 (by rwa [oc0] at hc)) hx,
        ← hL.r_next _ _ (by rwa [oc0] at hc) R1, ← md0, Nat.mod_add_mod, Nat.add_assoc]
    · rw [D.rot'_f1] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_gf _ _)
        (by simp only [fin3_add_two_add_two]; exact D.ne_gg (FlipAt.s1_ne _).symm)
      rw [(lσ x1 x2).2.2.2, rf1, D.oc_f1, W1, hrl, ← hL.r_next _ _ (by rwa [D.oc_f1] at hc) hx,
        Nat.mod_add_mod, Wl]
      congr 1; ring
    · rw [D.rot'_f2] at h; cases h
      rw [r2, oc2, W2, rg1, hrk]; congr 1; ring
    · rw [D.rot'_g0] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (D.ne_gf _ _) (D.ne_gg (FlipAt.s2_ne _).symm)
      rw [(lσ x1 x2).2.2.2, q0, od0, V0, ← hL.r_next _ _ (inp1 (by rwa [od0] at hc)) hx,
        ← hL.r_next _ _ (by rwa [od0] at hc) R2, ← md1, Nat.mod_add_mod, Nat.add_assoc]
    · rw [D.rot'_g1] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx
        (by simp only [fin3_add_two_add_two]; exact D.ne_ff (FlipAt.s1_ne _).symm) (D.ne_fg _ _)
      rw [(lσ x1 x2).2.2.2, rg1, D.oc_g1, V1, hrk, ← hL.r_next _ _ (by rwa [D.oc_g1] at hc) hx,
        Nat.mod_add_mod, Wk]
      congr 1; ring
    · rw [D.rot'_g2] at h; cases h
      rw [q2, od2, V2, rf1, hrl]; congr 1; ring
    · rw [D.rot'_O hf hg] at h
      obtain ⟨x, hx, rfl⟩ := Option.map_eq_some_iff.mp h
      obtain ⟨x1, x2⟩ := nold hx (fun e => hf (by rw [(Prod.mk.inj e).1]))
        (fun e => hg (by rw [(Prod.mk.inj e).1]))
      rw [(lσ x1 x2).2.2.2, rO c hf hg, oO c hf hg, WO c hf hg,
        hL.r_next _ _ (by rwa [oO c hf hg] at hc) hx]
  · -- r_last
    show (D.labels L).r c + (C.flip D).w c + 1 = L.deg (D.oc c)
    change L.input (D.oc c) at hc
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [D.rot'_f0, Option.map_eq_none_iff] at h
      rw [oc0] at hc ⊢
      have i0 := inp0 hc
      have hm := hL.md_bdry _ i0 h
      have e1 := hL.r_next _ _ hc R1
      rw [← md0, hm, Nat.mod_zero] at e1
      have e2 := hL.r_last _ i0 h
      rw [r0, W0, ← dg0]; omega
    · rw [D.rot'_f1, Option.map_eq_none_iff] at h
      rw [D.oc_f1] at hc ⊢
      have hm := hL.md_bdry _ hc h
      have e2 := hL.r_last _ hc h
      rw [rf1, hrl, hm, Nat.mod_zero, W1]; rw [Wl] at e2; omega
    · rw [D.rot'_f2] at h; cases h
    · rw [D.rot'_g0, Option.map_eq_none_iff] at h
      rw [od0] at hc ⊢
      have i0 := inp1 hc
      have hm := hL.md_bdry _ i0 h
      have e1 := hL.r_next _ _ hc R2
      rw [← md1, hm, Nat.mod_zero] at e1
      have e2 := hL.r_last _ i0 h
      rw [q0, V0, ← dg1]; omega
    · rw [D.rot'_g1, Option.map_eq_none_iff] at h
      rw [D.oc_g1] at hc ⊢
      have hm := hL.md_bdry _ hc h
      have e2 := hL.r_last _ hc h
      rw [rg1, hrk, hm, Nat.mod_zero, V1]; rw [Wk] at e2; omega
    · rw [D.rot'_g2] at h; cases h
    · rw [D.rot'_O hf hg, Option.map_eq_none_iff] at h
      rw [oO c hf hg] at hc ⊢
      rw [rO c hf hg, WO c hf hg]; exact hL.r_last _ hc h
  · -- r_first
    show (D.labels L).r c = 0
    change L.input (D.oc c) at hc
    rcases D.corner_cases c with rfl | rfl | rfl | rfl | rfl | rfl | ⟨hf, hg⟩
    · rw [D.twin'_f0, Option.map_eq_none_iff] at h
      rw [oc0] at hc; rw [r0]; exact hL.r_first _ hc h
    · rw [D.twin_f1] at h; cases h
    · rw [D.twin'_f2, Option.map_eq_none_iff] at h
      rw [oc2] at hc; rw [r2]; exact hL.r_first _ hc h
    · rw [D.twin'_g0, Option.map_eq_none_iff] at h
      rw [od0] at hc; rw [q0]; exact hL.r_first _ hc h
    · rw [D.twin_g1] at h; cases h
    · rw [D.twin'_g2, Option.map_eq_none_iff] at h
      rw [od2] at hc; rw [q2]; exact hL.r_first _ hc h
    · rw [D.twin'_O hf hg, Option.map_eq_none_iff] at h
      rw [oO c hf hg] at hc; rw [rO c hf hg]; exact hL.r_first _ hc h

end Cx

end DDG.NormalCoords
