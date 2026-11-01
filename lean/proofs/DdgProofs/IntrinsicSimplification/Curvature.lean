import DdgProofs.IntrinsicSimplification.Flatten

/-!
# Intrinsic simplification: flattening conserves curvature; the redistribution weights (Eq. 9)

`crates/ddg-intrinsic-simplification/SPEC.md` §3.4 (Commit), §4.2, S11, S12, §13; the paper
Sec. 4.1, Sec. 5.4, Eq. 9, Fig. 10.

**The mesh.** A Δ-complex is modelled by its faces `F` (finite) with the vertex labels of their
corners `corner f k ∈ V` (repetitions allowed) and their side lengths `len f k > 0` (side `k`
opposite corner `k`). The angle sum `Θ_v` adds the (extended) angles of every corner labelled
`v`. Flattening `i` with the scale `u` multiplies side `k` of a face by `e^{u m/2}`, `m` the
number of its two ends labelled `i` (Eq. 3 with `u_j = 0` for `j ≠ i`).

* **`iem-flatten-conserves-curvature`** — only `Θ` of `i` and of the vertices of its star change
  (`flatten_angleSum_eq_of_not_near`); over any set `T` of vertices containing every corner of the
  star, the total angle is unchanged (`flatten_sum_angleSum`; each face's angles sum to `π` before
  and after, extended angles included); hence with `K_v = c_v − Θ_v` (`c_v = 2π` inside, `π` on
  the boundary, mixed freely) and `i` flattened (`Θ̃_i = c_i`):
  `Σ_{j ∈ N_i} K̃_j = Σ_{j ∈ N_i} K_j + K_i` (`flatten_curvature_conserved`).
* **`iem-redistribution-weights`** — Eq. 9 with the flat-vertex rule of S11 (`weights`): the
  weights are non-negative and sum to one, and for a vertex that is not flat
  (`|K_i| > tol ≥ 0`) the denominator `Σ_l |K̃_l − K_l|` is at least `|K_i| > 0`
  (`weights_nonneg`, `weights_sum`, `abs_le_sum_abs_dK`).
-/

noncomputable section

namespace DdgProofs

namespace IntrinsicSimplification

open Real Finset

/-! ## `iem-flatten-conserves-curvature` -/

variable {V F : Type*} [DecidableEq V] [Fintype F]

/-- A Δ-complex given by the corners and side lengths of its faces. -/
structure FaceMesh (V F : Type*) where
  /-- vertex label of corner `k` of face `f` -/
  corner : F → Fin 3 → V
  /-- side `k` of face `f`, opposite corner `k` -/
  len : F → Fin 3 → ℝ

namespace FaceMesh

variable (M : FaceMesh V F)

/-- The angle at corner `k` of a face with sides `l` (extended angle). -/
def angle (l : Fin 3 → ℝ) (k : Fin 3) : ℝ := cornerAngle (l k) (l (k + 1)) (l (k + 2))

/-- The angle sum at `v`: all corners labelled `v`. -/
def angleSum (v : V) : ℝ := ∑ f, ∑ k, if M.corner f k = v then angle (M.len f) k else 0

/-- `1` if corner `k` of `f` is `i`. -/
def mult (i : V) (f : F) (k : Fin 3) : ℝ := if M.corner f k = i then 1 else 0

/-- Flattening `i` with scale `u` (Eq. 3, `u_j = 0` for `j ≠ i`). -/
def flatten (i : V) (u : ℝ) : FaceMesh V F where
  corner := M.corner
  len := fun f k => M.len f k * exp (u * (M.mult i f (k + 1) + M.mult i f (k + 2)) / 2)

/-- `f` is a face of the star of `i`. -/
def InStar (i : V) (f : F) : Prop := ∃ k, M.corner f k = i

/-- The angles of a face with positive sides sum to `π` (extended angles). -/
theorem sum_angle {l : Fin 3 → ℝ} (hl : ∀ k, 0 < l k) : ∑ k, angle l k = π := by
  rw [Fin.sum_univ_three]
  exact cornerAngle_sum_pos (hl 0) (hl 1) (hl 2)

omit [Fintype F] in
theorem flatten_len_pos (hl : ∀ f k, 0 < M.len f k) (i : V) (u : ℝ) :
    ∀ f k, 0 < (M.flatten i u).len f k :=
  fun f k => mul_pos (hl f k) (exp_pos _)

omit [Fintype F] in
/-- A face without a corner at `i` keeps its sides. -/
theorem flatten_len_of_not_inStar {i : V} {f : F} (h : ¬ M.InStar i f) (u : ℝ) :
    (M.flatten i u).len f = M.len f := by
  funext k
  have hm : ∀ k, M.mult i f k = 0 := fun k => by
    simp only [mult]; split_ifs with hk
    · exact absurd ⟨k, hk⟩ h
    · rfl
  simp [flatten, hm]

/-- **Only the star changes**: a vertex that is not a corner of any face of the star of `i` keeps
its angle sum.

[formula: iem-flatten-conserves-curvature] -/
theorem flatten_angleSum_eq_of_not_near {i v : V}
    (hv : ∀ f, M.InStar i f → ∀ k, M.corner f k ≠ v) (u : ℝ) :
    (M.flatten i u).angleSum v = M.angleSum v := by
  unfold angleSum
  apply Finset.sum_congr rfl
  intro f _
  by_cases hf : M.InStar i f
  · apply Finset.sum_congr rfl
    intro k _
    have : M.corner f k ≠ v := hv f hf k
    simp [flatten, this]
  · rw [M.flatten_len_of_not_inStar hf]
    rfl

/-- **The total angle over the star is unchanged**: for any set `T` of vertices that contains every
corner of every face of the star of `i`, `Σ_{v ∈ T} Θ̃_v = Σ_{v ∈ T} Θ_v` (each face's angles sum
to `π` before and after).

[formula: iem-flatten-conserves-curvature] -/
theorem flatten_sum_angleSum (hl : ∀ f k, 0 < M.len f k) {i : V} (T : Finset V)
    (hT : ∀ f, M.InStar i f → ∀ k, M.corner f k ∈ T) (u : ℝ) :
    ∑ v ∈ T, (M.flatten i u).angleSum v = ∑ v ∈ T, M.angleSum v := by
  -- exchange the sums: each face contributes the angles of its corners labelled in `T`
  have swap : ∀ N : FaceMesh V F, N.corner = M.corner → ∑ v ∈ T, N.angleSum v =
      ∑ f, ∑ k, if M.corner f k ∈ T then angle (N.len f) k else 0 := by
    intro N hN
    unfold angleSum
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl; intro f _
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl; intro k _
    rw [hN, Finset.sum_ite_eq]
  rw [swap (M.flatten i u) rfl, swap M rfl]
  apply Finset.sum_congr rfl
  intro f _
  by_cases hf : M.InStar i f
  · have hall : ∀ k, M.corner f k ∈ T := hT f hf
    simp only [hall, ite_true]
    rw [sum_angle (M.flatten_len_pos hl i u f), sum_angle (hl f)]
  · rw [M.flatten_len_of_not_inStar hf]

/-- **`iem-flatten-conserves-curvature`** (Fig. 10, local Gauss–Bonnet). Let `K_v = c_v − Θ_v`
(`c_v = 2π` at an interior vertex, `π` at a boundary vertex: Eq. 1 and Eq. 2 may be mixed), `N`
the neighbours of `i` (a set not containing `i` that, with `i`, contains every corner of the star)
and `u` a flattening scale with `Θ̃_i = c_i`. Then `Σ_{j ∈ N} K̃_j = Σ_{j ∈ N} K_j + K_i`.

[formula: iem-flatten-conserves-curvature] -/
theorem flatten_curvature_conserved (hl : ∀ f k, 0 < M.len f k) (c : V → ℝ) {i : V}
    (N : Finset V) (hi : i ∉ N) (hN : ∀ f, M.InStar i f → ∀ k, M.corner f k ∈ insert i N)
    {u : ℝ} (hflat : (M.flatten i u).angleSum i = c i) :
    ∑ j ∈ N, (c j - (M.flatten i u).angleSum j) =
      ∑ j ∈ N, (c j - M.angleSum j) + (c i - M.angleSum i) := by
  have h := M.flatten_sum_angleSum hl (insert i N) hN u
  rw [Finset.sum_insert hi, Finset.sum_insert hi, hflat] at h
  simp only [Finset.sum_sub_distrib]
  linarith

end FaceMesh

/-! ## `iem-redistribution-weights` (Eq. 9) -/

variable {ι : Type*} [Fintype ι]

/-- Eq. 9 with the flat-vertex rule of S11: `dK l = K̃_l − K_l` (or `κ̃ − κ` at a boundary
neighbour) over the distinct neighbours `ι = N_i`; if `|K_i| ≤ tol` the vertex is already flat
and every neighbour gets `1/|N_i|`. -/
def weights (dK : ι → ℝ) (Ki tol : ℝ) (j : ι) : ℝ :=
  if |Ki| ≤ tol then 1 / Fintype.card ι else |dK j| / ∑ l, |dK l|

/-- **The denominator of Eq. 9 is at least `|K_i|`**: by conservation (`Σ dK = K_i`,
`iem-flatten-conserves-curvature`) and the triangle inequality.

[formula: iem-redistribution-weights] -/
theorem abs_le_sum_abs_dK {dK : ι → ℝ} {Ki : ℝ} (hcons : ∑ l, dK l = Ki) :
    |Ki| ≤ ∑ l, |dK l| := by
  rw [← hcons]; exact Finset.abs_sum_le_sum_abs _ _

/-- The weights are non-negative.

[formula: iem-redistribution-weights] -/
theorem weights_nonneg (dK : ι → ℝ) (Ki tol : ℝ) (j : ι) : 0 ≤ weights dK Ki tol j := by
  unfold weights
  split_ifs
  · positivity
  · exact div_nonneg (abs_nonneg _) (Finset.sum_nonneg fun _ _ => abs_nonneg _)

/-- **The weights sum to one** (given conservation, `0 ≤ tol`, and at least one neighbour).

[formula: iem-redistribution-weights] -/
theorem weights_sum [Nonempty ι] {dK : ι → ℝ} {Ki tol : ℝ} (htol : 0 ≤ tol)
    (hcons : ∑ l, dK l = Ki) : ∑ j, weights dK Ki tol j = 1 := by
  unfold weights
  split_ifs with h
  · rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    have : (0 : ℝ) < Fintype.card ι := by exact_mod_cast Fintype.card_pos
    field_simp
  · have hpos : 0 < ∑ l, |dK l| := by
      have := abs_le_sum_abs_dK hcons
      push Not at h
      linarith
    rw [← Finset.sum_div, div_self hpos.ne']

/-- For a vertex that is not flat the weights are Eq. 9 itself, with a positive denominator.

[formula: iem-redistribution-weights] -/
theorem weights_eq_eq9 {dK : ι → ℝ} {Ki tol : ℝ} (htol : 0 ≤ tol) (hcons : ∑ l, dK l = Ki)
    (h : tol < |Ki|) (j : ι) :
    weights dK Ki tol j = |dK j| / ∑ l, |dK l| ∧ 0 < ∑ l, |dK l| := by
  refine ⟨by simp [weights, not_le.mpr h], ?_⟩
  have := abs_le_sum_abs_dK hcons
  linarith

end IntrinsicSimplification

end DdgProofs
