import DdgProofs.Developability.Basic

/-!
# Developability: midpoint 4-1 subdivision preserves the energies

`crates/ddg-developability/SPEC.md` §5.2, D20 (the paper's inset of Sec. 4.2: regular 4-1 split,
no smoothing). A face `t = (a, b, c)` with the edge midpoints `m_ab, m_bc, m_ca`
(`Ddg.Dev.midpoint`) becomes `(a, m_ab, m_ca)`, `(b, m_bc, m_ab)`, `(c, m_ca, m_bc)`,
`(m_ab, m_bc, m_ca)`. Proved:

* every child has the normal vector `n/4`, hence the unit normal of its parent (`child_normals`,
  `child_unit_normals`);
* the corner of a child at an old vertex has the angle of the parent's corner
  (`child_angles`), and the vector area of an old star is divided by `4`, so its vertex normal is
  kept (`vectorArea_children`, `vertexNormal_children`). Since all six energies of a star are
  functions of the normals and angles of its faces and of its vertex normal, the energy of every
  old vertex is unchanged;
* a new vertex (an edge midpoint) is surrounded by children of the one or two faces of its edge,
  so its normals take at most two values `N_1, N_2` (two planar arcs); for any star whose normals
  are all orthogonal to one unit vector `w` (e.g. `w = N_1 × N_2` normalized), the objectives of
  Eq. 3 and Eq. 6 vanish at `w`, the vertex normal and every log-mapped normal are orthogonal to
  `w`, so the intrinsic objectives vanish at `w` too (`coplanar_energies_vanish`); and a part with
  constant normal has width `0` (`partWidth_const`). With the energies being minima of
  nonnegative objectives (`dev-covariance-matrix`, `dev-intrinsic-width`), the energy of a new
  vertex is `0`: a hinge.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open scoped DdgProofs.Oracle
open Real Ddg.Dev

theorem midpoint_coords (a b : V3) :
    midpoint a b = ⟨(a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2⟩ := by
  unfold Ddg.Dev.midpoint
  ext <;> simp [Oracle.one_eq, Oracle.two_eq] <;> ring

/-! ## Children of a face -/

/-- **Children are coplanar with the parent**: each of the four children of `(a, b, c)` has the
normal vector `n/4`.

[formula: dev-midpoint-subdivision] -/
theorem child_normals (a b c : V3) :
    faceNormal a (midpoint a b) (midpoint c a) = smul (1 / 4) (faceNormal a b c) ∧
    faceNormal b (midpoint b c) (midpoint a b) = smul (1 / 4) (faceNormal a b c) ∧
    faceNormal c (midpoint c a) (midpoint b c) = smul (1 / 4) (faceNormal a b c) ∧
    faceNormal (midpoint a b) (midpoint b c) (midpoint c a) = smul (1 / 4) (faceNormal a b c) := by
  simp only [midpoint_coords]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> ext <;> simp [faceNormal_def] <;> ring

theorem normalize_smul_pos {s : ℝ} (hs : 0 < s) (n : V3) (hn : n ≠ zero) :
    normalize (smul s n) = normalize n := by
  rw [normalize_def, normalize_def, norm_smul_of_nonneg n hs.le]
  have hp := (norm_pos_iff' n).mpr hn
  ext <;> simp only [smul_x, smul_y, smul_z] <;> field_simp

/-- The children have the unit normal of the parent.

[formula: dev-midpoint-subdivision] -/
theorem child_unit_normals (a b c : V3) (h : faceNormal a b c ≠ zero) :
    faceUnitNormal a (midpoint a b) (midpoint c a) = faceUnitNormal a b c ∧
    faceUnitNormal b (midpoint b c) (midpoint a b) = faceUnitNormal a b c ∧
    faceUnitNormal c (midpoint c a) (midpoint b c) = faceUnitNormal a b c ∧
    faceUnitNormal (midpoint a b) (midpoint b c) (midpoint c a) = faceUnitNormal a b c := by
  obtain ⟨h1, h2, h3, h4⟩ := child_normals a b c
  simp only [faceUnitNormal_def, h1, h2, h3, h4, normalize_smul_pos (by norm_num : (0 : ℝ) < 1 / 4) _ h,
    and_self]

/-- `cornerAngleAt` only depends on the two edge vectors, and is invariant under scaling them by
the same positive factor. -/
theorem cornerAngleAt_scale {s : ℝ} (hs : 0 < s) (p q r p' q' r' : V3)
    (h1 : sub q' p' = smul s (sub q p)) (h2 : sub r' p' = smul s (sub r p)) :
    cornerAngleAt p' q' r' = cornerAngleAt p q r := by
  unfold cornerAngleAt
  rw [h1, h2]
  have ec : cross (smul s (sub q p)) (smul s (sub r p)) = smul (s * s) (cross (sub q p) (sub r p)) := by
    ext <;> simp <;> ring
  have ed : dot (smul s (sub q p)) (smul s (sub r p)) = (s * s) * dot (sub q p) (sub r p) := by
    simp only [dot_smul_l, dot_smul_r]; ring
  rw [ec, ed, norm_smul_of_nonneg _ (mul_pos hs hs).le]
  exact atan2_smul (mul_pos hs hs) _ _

/-- **Old angles are kept**: the corner of a child at an old vertex has the angle of the parent's
corner there.

[formula: dev-midpoint-subdivision] -/
theorem child_angles (a b c : V3) :
    cornerAngleAt a (midpoint a b) (midpoint c a) = cornerAngleAt a b c ∧
    cornerAngleAt b (midpoint b c) (midpoint a b) = cornerAngleAt b c a ∧
    cornerAngleAt c (midpoint c a) (midpoint b c) = cornerAngleAt c a b := by
  refine ⟨?_, ?_, ?_⟩ <;> apply cornerAngleAt_scale (by norm_num : (0 : ℝ) < 1 / 2) <;>
    simp only [midpoint_coords] <;> ext <;> simp <;> ring

/-- **Vertex normals of old vertices are kept**: replacing every face `(a_m, b_m, c_m)` of a star
(at the vertex `a_m`) by its child `(a_m, m_ab, m_ca)` divides the vector area by `4`.

[formula: dev-midpoint-subdivision] -/
theorem vectorArea_children (k : ℕ) (a b c : ℕ → V3) :
    vectorArea k (fun m => (a m, midpoint (a m) (b m), midpoint (c m) (a m))) =
      smul (1 / 4) (vectorArea k (fun m => (a m, b m, c m))) := by
  have h := fun m => (child_normals (a m) (b m) (c m)).1
  unfold vectorArea vsumTo
  simp only [h]
  ext <;> simp only [smul_x, smul_y, smul_z, sumTo_eq_sum_range, ← Finset.mul_sum] <;> ring

theorem vertexNormal_children (k : ℕ) (a b c : ℕ → V3)
    (h : vectorArea k (fun m => (a m, b m, c m)) ≠ zero) :
    vertexNormal k (fun m => (a m, midpoint (a m) (b m), midpoint (c m) (a m))) =
      vertexNormal k (fun m => (a m, b m, c m)) := by
  unfold vertexNormal
  rw [vectorArea_children, normalize_smul_pos (by norm_num) _ h]

/-! ## Stars with coplanar normals (the new stars) -/

theorem sumTo_eq_zero {n : ℕ} {f : ℕ → ℝ} (h : ∀ m, m < n → f m = 0) : Ddg.sumTo n f = 0 := by
  rw [sumTo_eq_sum_range]
  exact Finset.sum_eq_zero fun m hm => h m (Finset.mem_range.mp hm)

theorem maxTo1_eq_zero {n : ℕ} {f : ℕ → ℝ} (h : ∀ m, m ≤ n → f m = 0) : maxTo1 n f = 0 := by
  induction n with
  | zero => exact h 0 le_rfl
  | succ n ih =>
    show max (maxTo1 n f) (f (n + 1)) = 0
    rw [ih (fun m hm => h m (by omega)), h (n + 1) le_rfl, max_self]

theorem minTo1_eq_zero {n : ℕ} {f : ℕ → ℝ} (h : ∀ m, m ≤ n → f m = 0) : minTo1 n f = 0 := by
  induction n with
  | zero => exact h 0 le_rfl
  | succ n ih =>
    show min (minTo1 n f) (f (n + 1)) = 0
    rw [ih (fun m hm => h m (by omega)), h (n + 1) le_rfl, min_self]

/-- **Coplanar normals give zero energies** (a hinge): if all normals of a star are orthogonal to
`w`, then the objective of Eq. 3 (`covObjective`) and of Eq. 6 (`maxCovObjective`) vanish at `w`;
if moreover all face normal vectors are, the vector area and the vertex normal `n0` are orthogonal
to `w`, every log-mapped normal `Ñ = log_{n0} N` is orthogonal to `w`, and the intrinsic objectives
(`covObjective` of the `Ñ`, `widthAlong`) vanish at `w` (a tangent vector at `n0`).

[formula: dev-midpoint-subdivision] -/
theorem coplanar_energies_vanish (k : ℕ) (θ : ℕ → ℝ) (N : ℕ → V3) (face : ℕ → V3 × V3 × V3)
    (w : V3) (hN : ∀ m, dot w (N m) = 0)
    (hn : ∀ m, dot w (faceNormal (face m).1 (face m).2.1 (face m).2.2) = 0) :
    covObjective k θ N w = 0 ∧ maxCovObjective k N w = 0 ∧
    dot w (vectorArea k face) = 0 ∧ dot w (vertexNormal k face) = 0 ∧
    (∀ m, dot w (logMap (vertexNormal k face) (N m)) = 0) ∧
    covObjective k θ (fun m => logMap (vertexNormal k face) (N m)) w = 0 ∧
    widthAlong k (fun m => logMap (vertexNormal k face) (N m)) w = 0 := by
  have hS : dot w (vectorArea k face) = 0 := by
    unfold vectorArea
    rw [dot_smul_r]
    have : dot w (vsumTo k fun m => faceNormal (face m).1 (face m).2.1 (face m).2.2) =
        Ddg.sumTo k (fun m => dot w (faceNormal (face m).1 (face m).2.1 (face m).2.2)) := by
      simp only [vsumTo, dot_def, sumTo_eq_sum_range, Finset.mul_sum, ← Finset.sum_add_distrib]
    rw [this, sumTo_eq_zero (fun m _ => hn m), mul_zero]
  have hv : dot w (vertexNormal k face) = 0 := by
    unfold vertexNormal; rw [normalize_def, dot_smul_r, hS, mul_zero]
  have hlog : ∀ m, dot w (logMap (vertexNormal k face) (N m)) = 0 := by
    intro m
    unfold logMap
    rw [dot_smul_r, dot_sub_r, dot_smul_r, hN m, hv]; ring
  refine ⟨?_, ?_, hS, hv, hlog, ?_, ?_⟩
  · exact sumTo_eq_zero fun m _ => by rw [hN m]; ring
  · exact maxTo1_eq_zero fun m _ => by rw [hN m]; ring
  · exact sumTo_eq_zero fun m _ => by rw [hlog m]; ring
  · unfold widthAlong
    rw [maxTo1_eq_zero (fun m _ => by rw [dot_comm, hlog m]),
      minTo1_eq_zero (fun m _ => by rw [dot_comm, hlog m]), sub_zero]

/-- A part of `n ≥ 1` faces with one normal has width `0` (Eq. 1): each arc of a hinge.

[formula: dev-midpoint-subdivision] -/
theorem partWidth_const {n : ℕ} (hn : 0 < n) (N : V3) : partWidth n (fun _ => N) = 0 := by
  have hmean : partMean n (fun _ => N) = N := by
    unfold partMean vsumTo
    have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
    ext <;> simp only [smul_x, smul_y, smul_z, sumTo_eq_sum_range, Finset.sum_const,
      Finset.card_range, nsmul_eq_mul, scalar_ofNat, Oracle.one_eq] <;> field_simp
  unfold partWidth
  rw [hmean, sumTo_eq_zero (fun m _ => by simp [norm2_def]), zero_div]

end Dev

end DdgProofs
