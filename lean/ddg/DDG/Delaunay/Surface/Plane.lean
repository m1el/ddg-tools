import DDG.Delaunay.Planar.Basic
import DDG.Intrinsic.Diamond

/-!
# Delaunay layer 4, planar toolkit: rigid motions and triangles

Points are pairs `ℝ × ℝ` with the Euclidean distance `dist2` and the signed double area `sda`
(`DDG/Geometry/Plane.lean`); hull membership is by the signs of `sda`
(`DDG/Delaunay/Planar/Basic.lean`). This file packages

* **rigid motions** `Rigid` (a rotation `(cs, sn)` followed by a translation, the `rigid` of
  `DDG/Intrinsic/Diamond.lean`) with composition, inverse, the isometry and orientation
  facts, and `rigidOfSegment`: the unique rigid motion taking a segment onto a congruent one
  (`rigidOfSegment_fst`, `rigidOfSegment_snd`, `Rigid.ext_of_two`);
* **triangles** `Tri := Fin 3 → ℝ × ℝ`, positively oriented (`Tri.Pos`), with their closed hull,
  interior and open sides, and the barycentric numerators `bary` (`sda` of the point with the two
  other corners): `mem_hull_iff`, `mem_int_iff`, `mem_side_iff`; a point of an open side is
  neither interior nor a corner; the image of a triangle under a rigid motion.

The lemmas about straight segments through triangles are in `Segment.lean`, the angle lemmas in
`Angle.lean`.
-/

noncomputable section

namespace DDG.Delaunay.Surface

open DDG.Intrinsic (rigid dist2_rigid sda_rigid)
open DDG.Delaunay.Planar
open Set

/-- A point of the plane. -/
abbrev Pt := ℝ × ℝ

/-! ## Rigid motions -/

/-- A rigid motion: the rotation with cosine `cs` and sine `sn`, followed by the translation
by `t` (`DDG.Intrinsic.rigid`). -/
@[ext]
structure Rigid where
  cs : ℝ
  sn : ℝ
  t : Pt
  norm : cs ^ 2 + sn ^ 2 = 1

namespace Rigid

instance : CoeFun Rigid (fun _ => Pt → Pt) := ⟨fun g => rigid g.cs g.sn g.t⟩

theorem apply_def (g : Rigid) (p : Pt) : g p = rigid g.cs g.sn g.t p := rfl

theorem apply_eq (g : Rigid) (p : Pt) :
    g p = (g.t.1 + g.cs * p.1 - g.sn * p.2, g.t.2 + g.sn * p.1 + g.cs * p.2) := rfl

theorem apply_fst (g : Rigid) (p : Pt) : (g p).1 = g.t.1 + g.cs * p.1 - g.sn * p.2 := rfl
theorem apply_snd (g : Rigid) (p : Pt) : (g p).2 = g.t.2 + g.sn * p.1 + g.cs * p.2 := rfl

theorem rigid_fst (cs sn : ℝ) (t p : Pt) : (rigid cs sn t p).1 = t.1 + cs * p.1 - sn * p.2 := rfl
theorem rigid_snd (cs sn : ℝ) (t p : Pt) : (rigid cs sn t p).2 = t.2 + sn * p.1 + cs * p.2 := rfl

/-- The identity. -/
def id : Rigid := ⟨1, 0, (0, 0), by norm_num⟩

@[simp] theorem id_apply (p : Pt) : id p = p := by
  ext
  · simp [rigid_fst, id]
  · simp [rigid_snd, id]

/-- Composition `g ∘ h` (first `h`, then `g`). -/
def comp (g h : Rigid) : Rigid :=
  ⟨g.cs * h.cs - g.sn * h.sn, g.sn * h.cs + g.cs * h.sn, g h.t, by
    have hg := g.norm
    have hh := h.norm
    linear_combination (h.cs ^ 2 + h.sn ^ 2) * hg + hh⟩

@[simp] theorem comp_apply (g h : Rigid) (p : Pt) : comp g h p = g (h p) := by
  ext
  · simp only [rigid_fst, rigid_snd, comp]; ring
  · simp only [rigid_fst, rigid_snd, comp]; ring

theorem comp_assoc (g h k : Rigid) : comp (comp g h) k = comp g (comp h k) := by
  ext
  · simp only [comp]; ring
  · simp only [comp]; ring
  · simp only [comp, rigid_fst, rigid_snd]; ring
  · simp only [comp, rigid_fst, rigid_snd]; ring

/-- The inverse motion. -/
def inv (g : Rigid) : Rigid :=
  ⟨g.cs, -g.sn, (-(g.cs * g.t.1) - g.sn * g.t.2, g.sn * g.t.1 - g.cs * g.t.2), by
    have := g.norm; linear_combination this⟩

@[simp] theorem inv_apply_apply (g : Rigid) (p : Pt) : g.inv (g p) = p := by
  have h := g.norm
  ext
  · simp only [rigid_fst, rigid_snd, inv]; linear_combination p.1 * h
  · simp only [rigid_fst, rigid_snd, inv]; linear_combination p.2 * h

@[simp] theorem apply_inv_apply (g : Rigid) (p : Pt) : g (g.inv p) = p := by
  have h := g.norm
  ext
  · simp only [rigid_fst, rigid_snd, inv]; linear_combination (p.1 - g.t.1) * h
  · simp only [rigid_fst, rigid_snd, inv]; linear_combination (p.2 - g.t.2) * h

theorem comp_id (g : Rigid) : g.comp id = g := by
  ext
  · simp [comp, id]
  · simp [comp, id]
  · simp [comp, id, rigid_fst]
  · simp [comp, id, rigid_snd]

theorem id_comp (g : Rigid) : id.comp g = g := by
  ext
  · simp [comp, id]
  · simp [comp, id]
  · simp [comp, id, rigid_fst]
  · simp [comp, id, rigid_snd]

theorem injective (g : Rigid) : Function.Injective (g : Pt → Pt) :=
  Function.LeftInverse.injective (inv_apply_apply g)

theorem surjective (g : Rigid) : Function.Surjective (g : Pt → Pt) :=
  Function.RightInverse.surjective (apply_inv_apply g)

@[simp] theorem apply_inj (g : Rigid) {p q : Pt} : g p = g q ↔ p = q := g.injective.eq_iff

theorem dist2_apply (g : Rigid) (p q : Pt) : dist2 (g p) (g q) = dist2 p q :=
  dist2_rigid g.norm g.t p q

theorem sda_apply (g : Rigid) (p q r : Pt) : sda (g p) (g q) (g r) = sda p q r :=
  sda_rigid g.norm g.t p q r

/-- A rigid motion is affine. -/
theorem apply_combo (g : Rigid) (p q : Pt) {a b : ℝ} (h : a + b = 1) :
    g (a • p + b • q) = a • g p + b • g q := by
  ext
  · simp only [rigid_fst, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    linear_combination (-g.t.1) * h
  · simp only [rigid_snd, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    linear_combination (-g.t.2) * h

theorem apply_add_smul (g : Rigid) (p w : Pt) (s : ℝ) :
    g (p + s • w) = g p + s • (g (p + w) - g p) := by
  ext
  · simp only [rigid_fst, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_sub, Prod.snd_sub]; ring
  · simp only [rigid_snd, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_sub, Prod.snd_sub]; ring

theorem continuous (g : Rigid) : Continuous (g : Pt → Pt) := by
  show Continuous (rigid g.cs g.sn g.t)
  unfold rigid
  fun_prop

/-- Two rigid motions that agree on two distinct points are equal. -/
theorem ext_of_two {g g' : Rigid} {A B : Pt} (hAB : A ≠ B) (hA : g A = g' A) (hB : g B = g' B) :
    g = g' := by
  obtain ⟨A1, A2⟩ := A
  obtain ⟨B1, B2⟩ := B
  have a1 := congrArg Prod.fst hA
  have a2 := congrArg Prod.snd hA
  have b1 := congrArg Prod.fst hB
  have b2 := congrArg Prod.snd hB
  simp only [rigid_fst, rigid_snd] at a1 a2 b1 b2
  -- the rotation parts agree on the nonzero vector `B − A`
  have e1 : (g.cs - g'.cs) * (B1 - A1) - (g.sn - g'.sn) * (B2 - A2) = 0 := by linarith
  have e2 : (g.sn - g'.sn) * (B1 - A1) + (g.cs - g'.cs) * (B2 - A2) = 0 := by linarith
  have hv : (B1 - A1) ^ 2 + (B2 - A2) ^ 2 ≠ 0 := by
    intro h0
    apply hAB
    have h1 : B1 - A1 = 0 := by nlinarith [sq_nonneg (B1 - A1), sq_nonneg (B2 - A2)]
    have h2 : B2 - A2 = 0 := by nlinarith [sq_nonneg (B1 - A1), sq_nonneg (B2 - A2)]
    simp only [Prod.mk.injEq]
    constructor <;> linarith
  have hc : (g.cs - g'.cs) * ((B1 - A1) ^ 2 + (B2 - A2) ^ 2) = 0 := by
    linear_combination (B1 - A1) * e1 + (B2 - A2) * e2
  have hs : (g.sn - g'.sn) * ((B1 - A1) ^ 2 + (B2 - A2) ^ 2) = 0 := by
    linear_combination (B1 - A1) * e2 - (B2 - A2) * e1
  have hc' : g.cs = g'.cs := by
    have := (mul_eq_zero.1 hc).resolve_right hv; linarith
  have hs' : g.sn = g'.sn := by
    have := (mul_eq_zero.1 hs).resolve_right hv; linarith
  ext
  · exact hc'
  · exact hs'
  · rw [hc', hs'] at a1; linarith
  · rw [hc', hs'] at a2; linarith

theorem inv_comp (g : Rigid) : g.inv.comp g = id :=
  ext_of_two (A := (0, 0)) (B := (1, 0)) (by simp) (by simp) (by simp)

theorem comp_inv (g : Rigid) : g.comp g.inv = id :=
  ext_of_two (A := (0, 0)) (B := (1, 0)) (by simp) (by simp) (by simp)

end Rigid

/-! ### The rigid motion taking a segment onto a congruent one -/

theorem sqDist2_ne_zero {A B : Pt} (hAB : A ≠ B) : sqDist2 A B ≠ 0 := by
  intro h0
  apply hAB
  have h1 : (B.1 - A.1) ^ 2 + (B.2 - A.2) ^ 2 = 0 := h0
  have h2 : B.1 - A.1 = 0 := by nlinarith [sq_nonneg (B.1 - A.1), sq_nonneg (B.2 - A.2)]
  have h3 : B.2 - A.2 = 0 := by nlinarith [sq_nonneg (B.1 - A.1), sq_nonneg (B.2 - A.2)]
  ext <;> linarith

theorem sqDist2_pos {A B : Pt} (hAB : A ≠ B) : 0 < sqDist2 A B :=
  lt_of_le_of_ne (sqDist2_nonneg A B) (sqDist2_ne_zero hAB).symm

/-- Lagrange's identity. -/
theorem dot2_sq_add_cross2_sq (u v : Pt) :
    dot2 u v ^ 2 + cross2 u v ^ 2 = (u.1 ^ 2 + u.2 ^ 2) * (v.1 ^ 2 + v.2 ^ 2) := by
  simp only [dot2, cross2]; ring

/-- The cosine of the rigid motion of `rigidOfSegment`. -/
def segCs (A B A' B' : Pt) : ℝ := dot2 (B - A) (B' - A') / sqDist2 A B

/-- The sine of the rigid motion of `rigidOfSegment`. -/
def segSn (A B A' B' : Pt) : ℝ := cross2 (B - A) (B' - A') / sqDist2 A B

theorem segCs_sq_add_segSn_sq {A B A' B' : Pt} (h : dist2 A B = dist2 A' B') (hAB : A ≠ B) :
    segCs A B A' B' ^ 2 + segSn A B A' B' ^ 2 = 1 := by
  have hq : sqDist2 A B ≠ 0 := sqDist2_ne_zero hAB
  have hd : sqDist2 A B = sqDist2 A' B' := by rw [← dist2_sq, ← dist2_sq, h]
  have hl := dot2_sq_add_cross2_sq (B - A) (B' - A')
  simp only [Prod.fst_sub, Prod.snd_sub] at hl
  have e1 : (B.1 - A.1) ^ 2 + (B.2 - A.2) ^ 2 = sqDist2 A B := rfl
  have e2 : (B'.1 - A'.1) ^ 2 + (B'.2 - A'.2) ^ 2 = sqDist2 A' B' := rfl
  rw [e1, e2, ← hd] at hl
  unfold segCs segSn
  rw [div_pow, div_pow, ← add_div, hl, div_eq_one_iff_eq (pow_ne_zero 2 hq)]
  ring

/-- The rigid motion with `A ↦ A'` and `B ↦ B'`, for `|AB| = |A'B'|`, `A ≠ B` (the identity
otherwise). -/
def rigidOfSegment (A B A' B' : Pt) : Rigid :=
  if h : dist2 A B = dist2 A' B' ∧ A ≠ B then
    ⟨segCs A B A' B', segSn A B A' B',
      (A'.1 - (segCs A B A' B' * A.1 - segSn A B A' B' * A.2),
        A'.2 - (segSn A B A' B' * A.1 + segCs A B A' B' * A.2)),
      segCs_sq_add_segSn_sq h.1 h.2⟩
  else Rigid.id

theorem rigidOfSegment_fst {A B A' B' : Pt} (h : dist2 A B = dist2 A' B') (hAB : A ≠ B) :
    rigidOfSegment A B A' B' A = A' := by
  simp only [rigidOfSegment, h, hAB, ne_eq, not_false_eq_true, and_self, ↓reduceDIte]
  ext
  · simp only [Rigid.rigid_fst]; ring
  · simp only [Rigid.rigid_snd]; ring

theorem rigidOfSegment_snd {A B A' B' : Pt} (h : dist2 A B = dist2 A' B') (hAB : A ≠ B) :
    rigidOfSegment A B A' B' B = B' := by
  have hq : sqDist2 A B ≠ 0 := sqDist2_ne_zero hAB
  have hd : sqDist2 A B = sqDist2 A' B' := by rw [← dist2_sq, ← dist2_sq, h]
  simp only [rigidOfSegment, h, hAB, ne_eq, not_false_eq_true, and_self, ↓reduceDIte]
  ext
  · simp only [Rigid.rigid_fst, segCs, segSn]
    simp only [sqDist2, dot2, cross2, Prod.fst_sub, Prod.snd_sub] at hd hq ⊢
    field_simp
    ring
  · simp only [Rigid.rigid_snd, segCs, segSn]
    simp only [sqDist2, dot2, cross2, Prod.fst_sub, Prod.snd_sub] at hd hq ⊢
    field_simp
    ring

/-- `rigidOfSegment` is the only rigid motion with `A ↦ A'`, `B ↦ B'`. -/
theorem rigidOfSegment_unique {A B A' B' : Pt} (h : dist2 A B = dist2 A' B') (hAB : A ≠ B)
    {g : Rigid} (hA : g A = A') (hB : g B = B') : g = rigidOfSegment A B A' B' :=
  Rigid.ext_of_two hAB (by rw [hA, rigidOfSegment_fst h hAB])
    (by rw [hB, rigidOfSegment_snd h hAB])

theorem rigidOfSegment_self {A B : Pt} (hAB : A ≠ B) : rigidOfSegment A B A B = Rigid.id :=
  (rigidOfSegment_unique rfl hAB (Rigid.id_apply A) (Rigid.id_apply B)).symm

/-! ## Triangles -/

/-- A triangle: its three corners. Side `k` runs from corner `k` to corner `k + 1`, opposite
corner `k + 2` (the convention of `Halfedge`). -/
abbrev Tri := Fin 3 → Pt

namespace Tri

variable (T : Tri)

/-- Positively oriented and nondegenerate. -/
def Pos : Prop := 0 < sda (T 0) (T 1) (T 2)

/-- The closed triangle. -/
def hull : Set Pt := convexHull ℝ {T 0, T 1, T 2}

/-- The open triangle. -/
def int : Set Pt := interior T.hull

/-- The open side `k`, from corner `k` to corner `k + 1`. -/
def side (k : Fin 3) : Set Pt := openSegment ℝ (T k) (T (k + 1))

/-- The barycentric numerator of corner `k` at `x`: `sda` of the two other corners (in the
cyclic order `k + 1, k + 2`) with `x`. It is the signed double area at `x = T k` and `0` on the
opposite side `k + 1`. -/
def bary (k : Fin 3) (x : Pt) : ℝ := sda (T (k + 1)) (T (k + 2)) x

/-- The image under a rigid motion. -/
def map (g : Rigid) : Tri := fun k => g (T k)

variable {T}

/-- Membership in the closed triangle: all three barycentric numerators are `≥ 0`. -/
theorem mem_hull_iff (hT : T.Pos) {x : Pt} : x ∈ T.hull ↔ ∀ k, 0 ≤ T.bary k x := by
  have h := mem_hull3_iff (a := T 0) (b := T 1) (c := T 2) hT.ne' (x := x)
  unfold hull
  rw [h, Fin.forall_fin_succ, Fin.forall_fin_two]
  simp only [bary, Fin.isValue, Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Fin.reduceAdd]
  rw [sda_rotate x (T 1) (T 2), ← sda_rotate (T 2) (T 0) x]
  have hT' : 0 < sda (T 0) (T 1) (T 2) := hT
  constructor
  · rintro ⟨h1, h2, h3⟩
    exact ⟨(mul_nonneg_iff_of_pos_left hT').1 h1, (mul_nonneg_iff_of_pos_left hT').1 h2,
      (mul_nonneg_iff_of_pos_left hT').1 h3⟩
  · rintro ⟨h1, h2, h3⟩
    exact ⟨mul_nonneg hT'.le h1, mul_nonneg hT'.le h2, mul_nonneg hT'.le h3⟩

/-- Membership in the open triangle: all three barycentric numerators are `> 0`. -/
theorem mem_int_iff (hT : T.Pos) {x : Pt} : x ∈ T.int ↔ ∀ k, 0 < T.bary k x := by
  have hne : sda (T 0) (T 1) (T 2) ≠ 0 := hT.ne'
  have hT' : 0 < sda (T 0) (T 1) (T 2) := hT
  rw [Fin.forall_fin_succ, Fin.forall_fin_two]
  simp only [bary, Fin.isValue, Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Fin.reduceAdd]
  rw [sda_rotate x (T 1) (T 2), ← sda_rotate (T 2) (T 0) x]
  constructor
  · intro hx
    obtain ⟨h1, h2, h3⟩ := pos_of_mem_interior_hull3 hne hx
    exact ⟨(mul_pos_iff_of_pos_left hT').1 h1, (mul_pos_iff_of_pos_left hT').1 h2,
      (mul_pos_iff_of_pos_left hT').1 h3⟩
  · rintro ⟨h1, h2, h3⟩
    exact mem_interior_hull3 hne (mul_pos hT' h1) (mul_pos hT' h2) (mul_pos hT' h3)

theorem int_subset_hull : T.int ⊆ T.hull := interior_subset

theorem isClosed_hull : IsClosed T.hull := isClosed_hull3 _ _ _

theorem isOpen_int : IsOpen T.int := isOpen_interior

theorem convex_hull : Convex ℝ T.hull := convex_convexHull ℝ _

theorem corner_mem_hull (k : Fin 3) : T k ∈ T.hull := by
  apply subset_convexHull
  fin_cases k <;> simp

theorem hull_subset_closure_int (hT : T.Pos) : T.hull ⊆ closure T.int :=
  hull3_subset_closure_interior hT.ne'

theorem int_nonempty (hT : T.Pos) : T.int.Nonempty := interior_hull3_nonempty hT.ne'

theorem bary_affine (k : Fin 3) (p w : Pt) (s : ℝ) :
    T.bary k (p + s • w) = T.bary k p + s * cross2 (T (k + 2) - T (k + 1)) w := by
  simp only [bary, sda, cross2, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    smul_eq_mul, Prod.fst_sub, Prod.snd_sub]
  ring

theorem bary_combo (k : Fin 3) (p q : Pt) {a b : ℝ} (h : a + b = 1) :
    T.bary k (a • p + b • q) = a * T.bary k p + b * T.bary k q := by
  simp only [bary, sda, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  have hb : b = 1 - a := by linarith
  subst hb
  ring

theorem continuous_bary (k : Fin 3) : Continuous (T.bary k) := continuous_sda₃ _ _

/-- The three numerators add up to the signed double area. -/
theorem sum_bary (x : Pt) : T.bary 0 x + T.bary 1 x + T.bary 2 x = sda (T 0) (T 1) (T 2) := by
  simp only [bary, sda, Fin.isValue, Fin.reduceAdd]
  ring

theorem bary_self (k : Fin 3) : T.bary k (T k) = sda (T 0) (T 1) (T 2) := by
  fin_cases k <;> simp only [bary, sda, Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one,
    Fin.reduceFinMk] <;> ring

/-- The numerator of a corner vanishes at every other corner. -/
theorem bary_corner_ne {j k : Fin 3} (h : j ≠ k) : T.bary k (T j) = 0 := by
  fin_cases j <;> fin_cases k <;>
    simp only [bary, sda, Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk] <;>
    first | ring1 | exact absurd rfl h

theorem bary_add_one (k : Fin 3) : T.bary k (T (k + 1)) = 0 :=
  bary_corner_ne (by omega)

theorem bary_add_two (k : Fin 3) : T.bary k (T (k + 2)) = 0 :=
  bary_corner_ne (by omega)

/-- Cramer's rule: a point is the barycentric combination of the corners. -/
theorem cramer (hT : T.Pos) (x : Pt) :
    x = (T.bary 0 x / sda (T 0) (T 1) (T 2)) • T 0 + (T.bary 1 x / sda (T 0) (T 1) (T 2)) • T 1 +
      (T.bary 2 x / sda (T 0) (T 1) (T 2)) • T 2 := by
  have h := sda_cramer (a := T 0) (b := T 1) (c := T 2) hT.ne' x
  simp only [bary, Fin.isValue, Fin.reduceAdd]
  rw [sda_rotate x (T 1) (T 2), ← sda_rotate (T 2) (T 0) x]
  exact h

/-- Two vanishing numerators: the point is the third corner. -/
theorem eq_corner_of_bary (hT : T.Pos) {k : Fin 3} {x : Pt} (h1 : T.bary (k + 1) x = 0)
    (h2 : T.bary (k + 2) x = 0) : x = T k := by
  have hs := sum_bary (T := T) x
  have hc := cramer hT x
  have hne : sda (T 0) (T 1) (T 2) ≠ 0 := hT.ne'
  fin_cases k <;> simp only [Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one,
    Fin.reduceFinMk] at h1 h2 ⊢
  · have h0 : T.bary 0 x = sda (T 0) (T 1) (T 2) := by linarith
    rw [h0, h1, h2, div_self hne] at hc
    simpa using hc
  · have h0 : T.bary 1 x = sda (T 0) (T 1) (T 2) := by linarith
    rw [h0, h1, h2, div_self hne] at hc
    simpa using hc
  · have h0 : T.bary 2 x = sda (T 0) (T 1) (T 2) := by linarith
    rw [h0, h1, h2, div_self hne] at hc
    simpa using hc

/-- The open side `k`: the numerator of the opposite corner `k + 2` vanishes, the two others are
positive. -/
theorem mem_side_iff (hT : T.Pos) {k : Fin 3} {x : Pt} :
    x ∈ T.side k ↔ T.bary (k + 2) x = 0 ∧ 0 < T.bary k x ∧ 0 < T.bary (k + 1) x := by
  have hT' : 0 < sda (T 0) (T 1) (T 2) := hT
  have hne : sda (T 0) (T 1) (T 2) ≠ 0 := hT.ne'
  constructor
  · intro hx
    rw [side, openSegment_eq_image] at hx
    obtain ⟨θ, ⟨h0, h1⟩, rfl⟩ := hx
    rw [bary_combo _ _ _ (by ring), bary_combo _ _ _ (by ring), bary_combo _ _ _ (by ring),
      bary_self, bary_add_one, bary_corner_ne (by omega), bary_corner_ne (by omega),
      bary_corner_ne (by omega), bary_self]
    refine ⟨by ring, ?_, ?_⟩ <;> nlinarith
  · rintro ⟨h0, h1, h2⟩
    rw [side, openSegment_eq_image]
    refine ⟨T.bary (k + 1) x / sda (T 0) (T 1) (T 2), ⟨div_pos h2 hT', ?_⟩, ?_⟩
    · rw [div_lt_one hT']
      have hs := sum_bary (T := T) x
      fin_cases k <;> simp only [Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one,
        Fin.reduceFinMk] at h0 h1 h2 ⊢ <;> linarith
    · have hc := cramer hT x
      have hs := sum_bary (T := T) x
      have hk : (1 - T.bary (k + 1) x / sda (T 0) (T 1) (T 2)) =
          T.bary k x / sda (T 0) (T 1) (T 2) := by
        field_simp
        fin_cases k <;> simp only [Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one,
          Fin.reduceFinMk] at h0 h1 h2 ⊢ <;> linarith
      show (1 - T.bary (k + 1) x / sda (T 0) (T 1) (T 2)) • T k +
        (T.bary (k + 1) x / sda (T 0) (T 1) (T 2)) • T (k + 1) = x
      rw [hk]
      fin_cases k <;> simp only [Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one,
        Fin.reduceFinMk] at h0 h1 h2 hc ⊢
      · rw [h0] at hc; simpa using hc.symm
      · rw [h0] at hc; simpa [add_comm] using hc.symm
      · rw [h0] at hc; simpa [add_assoc, add_comm] using hc.symm

theorem side_subset_hull (hT : T.Pos) (k : Fin 3) : T.side k ⊆ T.hull := by
  intro x hx
  rw [mem_side_iff hT] at hx
  rw [mem_hull_iff hT]
  intro j
  rcases (by omega : j = k ∨ j = k + 1 ∨ j = k + 2) with rfl | rfl | rfl
  · exact hx.2.1.le
  · exact hx.2.2.le
  · exact le_of_eq hx.1.symm

theorem side_symm (T : Tri) (k : Fin 3) : T.side k = openSegment ℝ (T (k + 1)) (T k) :=
  openSegment_symm ℝ _ _

/-- A point of an open side is not interior. -/
theorem not_mem_int_of_mem_side (hT : T.Pos) {k : Fin 3} {x : Pt} (hx : x ∈ T.side k) :
    x ∉ T.int := by
  intro h
  rw [mem_int_iff hT] at h
  rw [mem_side_iff hT] at hx
  linarith [h (k + 2), hx.1]

/-- A point of an open side is not a corner. -/
theorem ne_corner_of_mem_side (hT : T.Pos) {k : Fin 3} {x : Pt} (hx : x ∈ T.side k) (j : Fin 3) :
    x ≠ T j := by
  rintro rfl
  rw [mem_side_iff hT] at hx
  obtain ⟨h0, h1, h2⟩ := hx
  have hT' : 0 < sda (T 0) (T 1) (T 2) := hT
  fin_cases k <;> fin_cases j <;>
    simp only [bary, Fin.isValue, Fin.reduceAdd, Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk]
      at h0 h1 h2 <;>
    simp only [sda] at h0 h1 h2 hT' <;> ring_nf at h0 h1 h2 hT' <;> linarith

/-- A point lies on at most one open side. -/
theorem side_eq_of_mem (hT : T.Pos) {k k' : Fin 3} {x : Pt} (hk : x ∈ T.side k)
    (hk' : x ∈ T.side k') : k = k' := by
  by_contra hne
  have h1 := (mem_side_iff hT).1 hk
  have h2 := (mem_side_iff hT).1 hk'
  -- two numerators vanish: `x` is the third corner
  rcases (by omega : k' = k + 1 ∨ k' = k + 2) with rfl | rfl
  · rw [show k + 1 + 2 = k by omega] at h2
    exact ne_corner_of_mem_side hT hk (k + 1)
      (eq_corner_of_bary hT (k := k + 1) (by rw [show k + 1 + 1 = k + 2 by omega]; exact h1.1)
        (by rw [show k + 1 + 2 = k by omega]; exact h2.1))
  · rw [show k + 2 + 2 = k + 1 by omega] at h2
    exact ne_corner_of_mem_side hT hk k
      (eq_corner_of_bary hT (k := k) h2.1 h1.1)

/-- A corner is not interior. -/
theorem corner_not_mem_int (hT : T.Pos) (k : Fin 3) : T k ∉ T.int := by
  intro h
  rw [mem_int_iff hT] at h
  have := h (k + 1)
  rw [bary_corner_ne (by omega)] at this
  exact lt_irrefl _ this

/-- A point of the closed triangle that is neither interior nor a corner lies on an open side. -/
theorem exists_side_of_mem_hull (hT : T.Pos) {x : Pt} (hx : x ∈ T.hull) (hi : x ∉ T.int)
    (hc : ∀ k, x ≠ T k) : ∃ k, x ∈ T.side k := by
  rw [mem_hull_iff hT] at hx
  rw [mem_int_iff hT] at hi
  push Not at hi
  obtain ⟨j, hj⟩ := hi
  have hj0 : T.bary j x = 0 := le_antisymm hj (hx j)
  refine ⟨j + 1, ?_⟩
  rw [mem_side_iff hT, show j + 1 + 2 = j by omega, show j + 1 + 1 = j + 2 by omega]
  refine ⟨hj0, ?_, ?_⟩
  · rcases (hx (j + 1)).lt_or_eq with h | h
    · exact h
    · exact absurd (eq_corner_of_bary hT (k := j + 2)
        (by rw [show j + 2 + 1 = j by omega]; exact hj0)
        (by rw [show j + 2 + 2 = j + 1 by omega]; exact h.symm)) (hc _)
  · rcases (hx (j + 2)).lt_or_eq with h | h
    · exact h
    · exact absurd (eq_corner_of_bary hT (k := j + 1)
        (by rw [show j + 1 + 1 = j + 2 by omega]; exact h.symm)
        (by rw [show j + 1 + 2 = j by omega]; exact hj0)) (hc _)

/-! ### Rigid images -/

theorem map_pos (g : Rigid) (hT : T.Pos) : (T.map g).Pos := by
  unfold Pos map
  rw [g.sda_apply]
  exact hT

theorem bary_map (g : Rigid) (k : Fin 3) (x : Pt) : (T.map g).bary k (g x) = T.bary k x := by
  unfold bary map
  exact g.sda_apply _ _ _

theorem mem_map_hull_iff (g : Rigid) (hT : T.Pos) {x : Pt} :
    g x ∈ (T.map g).hull ↔ x ∈ T.hull := by
  rw [mem_hull_iff (map_pos g hT), mem_hull_iff hT]
  simp only [bary_map]

theorem mem_map_int_iff (g : Rigid) (hT : T.Pos) {x : Pt} :
    g x ∈ (T.map g).int ↔ x ∈ T.int := by
  rw [mem_int_iff (map_pos g hT), mem_int_iff hT]
  simp only [bary_map]

theorem mem_map_side_iff (g : Rigid) (hT : T.Pos) {k : Fin 3} {x : Pt} :
    g x ∈ (T.map g).side k ↔ x ∈ T.side k := by
  rw [mem_side_iff (map_pos g hT), mem_side_iff hT]
  simp only [bary_map]

theorem map_hull (g : Rigid) (hT : T.Pos) : (T.map g).hull = g '' T.hull := by
  ext y
  obtain ⟨x, rfl⟩ := g.surjective y
  rw [mem_map_hull_iff g hT, g.injective.mem_set_image]

theorem map_int (g : Rigid) (hT : T.Pos) : (T.map g).int = g '' T.int := by
  ext y
  obtain ⟨x, rfl⟩ := g.surjective y
  rw [mem_map_int_iff g hT, g.injective.mem_set_image]

theorem map_side (g : Rigid) (hT : T.Pos) (k : Fin 3) : (T.map g).side k = g '' T.side k := by
  ext y
  obtain ⟨x, rfl⟩ := g.surjective y
  rw [mem_map_side_iff g hT, g.injective.mem_set_image]

theorem map_map (g h : Rigid) : (T.map h).map g = T.map (g.comp h) := by
  funext k; simp [map]

@[simp] theorem map_id : T.map Rigid.id = T := by funext k; simp [map]

end Tri

end DDG.Delaunay.Surface
