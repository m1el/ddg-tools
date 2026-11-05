import Mathlib

/-!
# Boundary First Flattening, step II: from boundary data to a closed plane curve

Statements about the formulas evaluated by `crates/ddg-bff/src/curve.rs` and by
`flatten_to_disk` / `hilbert_transform` in `crates/ddg-bff/src/disk.rs` (Sawhney & Crane,
*Boundary First Flattening*, ACM TOG 2017, Sec. 4.5 and 6.5, Eq. 17, 19, 20). Every Rust
expression is transcribed as a `def`; the theorems are about these definitions over `ℝ` and `ℂ`
(exact arithmetic, not floating point).

## Conventions

A boundary loop has vertices `0, …, n − 1` in order with the surface on the left; edge `i` runs
from vertex `i` to vertex `i + 1 (mod n)`; the exterior angle `k_i` at vertex `i` is the turn from
edge `i − 1` to edge `i`, positive to the left.

* Statements that only need the cyclic successor are made for a finite type `β` of boundary
  positions and a permutation `next : Equiv.Perm β` (`prev = next.symm`). The Rust indexing is
  `β = Fin n`, `next = finRotate n`: `finRotate_val` and `finRotate_symm_val` say that these are
  `i ↦ (i + 1) % n` and `j ↦ (j + n − 1) % n`. Several loops at once are a permutation with
  several cycles.
* The sequential constructions (tangents, positions) are indexed by `ℕ`, with sums over
  `Finset.range`. `Fin.sum_univ_eq_sum_range` converts between the two.
* Plane vectors are complex numbers: `(x, y) = x + i y`, a unit tangent is `exp (i φ)`.

## Formula ids and results

* `bff-target-boundary-length` (a definition; it is cited, no theorem here carries its tag):
  `targetBoundaryLength`, positive for positive lengths (`targetBoundaryLength_pos`), the identity
  for `u = 0` (`targetBoundaryLength_zero`).
* `bff-dual-boundary-length`: `dualBoundaryLength`; the dual lengths have the same sum as the
  edge lengths (`sum_dualBoundaryLength`) and are positive for positive lengths
  (`dualBoundaryLength_pos`).
* `bff-close-exterior-angles`: `closedExteriorAngle`; afterwards the angles sum to `2π`
  (`sum_closedExteriorAngle`); nothing changes if the defect is zero
  (`closedExteriorAngle_of_defect_eq_zero`).
* `bff-uniformization-angles`: `uniformizationAngle`; the angles sum to `2π`
  (`sum_uniformizationAngle`) and are positive (`uniformizationAngle_pos`); with the dual lengths
  of the target lengths, as in the code: `flattenToDisk_angles`.
* `bff-boundary-tangents`: `tangentAngle`, `boundaryTangent`; unit length
  (`norm_boundaryTangent`), the turn at vertex `i ≥ 1` is `k_i` (`boundaryTangent_succ`), the turn
  at vertex `0` is `k_0` exactly when `Σ k ∈ 2πℤ` (`boundaryTangent_closes_iff`,
  `boundaryTangent_closes`, `boundaryTangent_turn`).
* `bff-boundary-curve`: `boundaryCurve`; edge `p` is `ℓ_p T_p` (`boundaryCurve_succ_sub`) of
  length `|ℓ_p|` (`norm_boundaryCurve_succ_sub`); the polygon closes iff `Σ ℓ_i T_i = 0`
  (`boundaryCurve_closingEdge_iff`). The gap-spreading step is `spreadCurve`: it starts and ends
  at `0` (`spreadCurve_zero`, `spreadCurve_end`) and changes every edge by its share of the gap
  (`spreadCurve_succ_sub`, `spreadCurve_closingEdge`).
* `bff-closing-lengths`: `closingMatrix` (the code's `S`), `closedLengths` (the code's update);
  the adjusted lengths satisfy the closing constraints iff `S lam = T ℓ*`
  (`closedLengths_closes_iff`), and then they are the unique minimizer of the weighted distance to
  the target lengths among all lengths that close (`closedLengths_optimal`,
  `closedLengths_unique`); the formula of the paper for invertible `S`
  (`closedLengthsPaper_isMinimizer`); one loop without shared lengths
  (`boundaryCurve_loopClosedLength`, `loopClosedLength_optimal`).
* `bff-hilbert-transform`: `hilbertRhs`; it is the gradient of the signed area
  (`loopSignedArea_eq_sum_hilbertRhs_mul`), sums to zero (`sum_hilbertRhs`), and `L b = rhs` holds
  exactly for the minimizers `b` of the conformal energy (`hilbertTransform_minimizes`,
  `hilbertTransform_of_minimizes`, `hilbertTransform_minimizes_conformalEnergy_iff`;
  `hilbertVector_minimizes_conformalEnergy_iff` for the vector that the code returns and a
  Laplacian that is `Matrix.PosSemidef`).

## Not covered

The rank-deficient branch of `close_group` (`solve_semidefinite` with a pivot threshold, dropping
constraint directions until the lengths change by an acceptable amount) is not modelled: the
statements below are about an exact solution `lam` of `S lam = T ℓ*`. Such a solution exists for
all positive weights, also when `S` is singular (`exists_closing_multiplier`), and every exact
solution gives the same lengths (`closedLengths_eq_of_solutions`); whether the pivoted
factorization with its threshold finds one is not addressed. The assembly of the constraint
matrix from the `sides` table and the grouping of coupled loops are not modelled either: the
general statement takes an arbitrary matrix `T`.
-/

noncomputable section

namespace DDG

open Finset Matrix
open scoped Real
open Complex (I)

/-! ## The Rust indexing of a boundary loop -/

section Indexing

/-- On `Fin n` the permutation `finRotate n` is the cyclic successor of the Rust code,
`i ↦ (i + 1) % n`. -/
theorem finRotate_val {n : ℕ} (i : Fin n) : ((finRotate n i : Fin n) : ℕ) = (i + 1) % n := by
  cases n with
  | zero => exact i.elim0
  | succ n =>
    rw [coe_finRotate]
    split_ifs with h
    · subst h; simp
    · have : (i : ℕ) < n := Fin.val_lt_last h
      rw [Nat.mod_eq_of_lt (by omega)]

/-- The inverse of `finRotate n` is the cyclic predecessor of the Rust code,
`j ↦ (j + n − 1) % n`. -/
theorem finRotate_symm_val {n : ℕ} (j : Fin n) :
    (((finRotate n).symm j : Fin n) : ℕ) = (j + n - 1) % n := by
  have hn : 0 < n := j.pos
  have hlt : (j + n - 1) % n < n := Nat.mod_lt _ hn
  have key : (finRotate n).symm j = ⟨(j + n - 1) % n, hlt⟩ := by
    rw [Equiv.symm_apply_eq]
    apply Fin.ext
    rw [finRotate_val]
    show (j : ℕ) = ((j + n - 1) % n + 1) % n
    rcases Nat.eq_zero_or_pos (j : ℕ) with h0 | hpos
    · have e : ((j : ℕ) + n - 1) % n = n - 1 := by
        rw [h0, Nat.zero_add]
        exact Nat.mod_eq_of_lt (by omega)
      rw [e, Nat.sub_add_cancel hn, Nat.mod_self, h0]
    · have e : ((j : ℕ) + n - 1) % n = j - 1 := by
        have e' : (j : ℕ) + n - 1 = (j - 1) + n := by omega
        rw [e', Nat.add_mod_right]
        exact Nat.mod_eq_of_lt (by omega)
      rw [e, Nat.sub_add_cancel hpos, Nat.mod_eq_of_lt j.isLt]
  rw [key]

end Indexing

/-! ## Target lengths (`target_boundary_lengths`)

Formula id `bff-target-boundary-length`. This is a definition (paper Alg. 1, line 9); only two
sanity properties are proved. -/

section Target

variable {β : Type*}

/-- Target length of boundary edge `i` (from vertex `i` to `next i`):
`ℓ*_i = e^{(u_i + u_{next i}) / 2} ℓ_i`, with `u` the log scale factors at the boundary vertices
and `ℓ` the edge lengths of the input metric. Rust:
`(0.5 * (scale_factors[i] + scale_factors[(i + 1) % n])).exp() * lengths[i]`. -/
def targetBoundaryLength (next : Equiv.Perm β) (u ℓ : β → ℝ) (i : β) : ℝ :=
  Real.exp (1 / 2 * (u i + u (next i))) * ℓ i

/-- A target length is positive if the input length is. -/
theorem targetBoundaryLength_pos (next : Equiv.Perm β) (u : β → ℝ) {ℓ : β → ℝ} {i : β}
    (h : 0 < ℓ i) : 0 < targetBoundaryLength next u ℓ i :=
  mul_pos (Real.exp_pos _) h

/-- For `u = 0` the target lengths are the input lengths. -/
theorem targetBoundaryLength_zero (next : Equiv.Perm β) (ℓ : β → ℝ) :
    targetBoundaryLength next 0 ℓ = ℓ := by
  funext i
  simp [targetBoundaryLength]

/-- `targetBoundaryLength` in the Rust indexing. -/
theorem targetBoundaryLength_finRotate {n : ℕ} (u ℓ : ℕ → ℝ) (i : Fin n) :
    targetBoundaryLength (finRotate n) (fun j => u j) (fun j => ℓ j) i =
      Real.exp (1 / 2 * (u i + u ((i + 1) % n))) * ℓ i := by
  simp only [targetBoundaryLength, finRotate_val]

end Target

/-! ## Dual lengths (`dual_boundary_lengths`) -/

section Dual

variable {β : Type*}

/-- Dual length of boundary vertex `j`: half of the two boundary edges at `j`,
`½ (ℓ_{prev j} + ℓ_j)` (edge `prev j` ends at `j`, edge `j` starts there). Rust:
`0.5 * (lengths[(j + n - 1) % n] + lengths[j])`. -/
def dualBoundaryLength (next : Equiv.Perm β) (ℓ : β → ℝ) (j : β) : ℝ :=
  1 / 2 * (ℓ (next.symm j) + ℓ j)

/-- **The dual lengths of the boundary vertices have the same sum as the lengths of the boundary
edges**: every edge gives one half to each of its two ends.

[formula: bff-dual-boundary-length] -/
theorem sum_dualBoundaryLength [Fintype β] (next : Equiv.Perm β) (ℓ : β → ℝ) :
    ∑ j, dualBoundaryLength next ℓ j = ∑ i, ℓ i := by
  simp only [dualBoundaryLength]
  rw [← Finset.mul_sum, Finset.sum_add_distrib, Equiv.sum_comp next.symm ℓ]
  ring

/-- The dual lengths are positive if the edge lengths are.

[formula: bff-dual-boundary-length] -/
theorem dualBoundaryLength_pos (next : Equiv.Perm β) {ℓ : β → ℝ} (hℓ : ∀ i, 0 < ℓ i) (j : β) :
    0 < dualBoundaryLength next ℓ j := by
  have h1 := hℓ (next.symm j)
  have h2 := hℓ j
  unfold dualBoundaryLength
  positivity

/-- `dualBoundaryLength` in the Rust indexing. -/
theorem dualBoundaryLength_finRotate {n : ℕ} (ℓ : ℕ → ℝ) (j : Fin n) :
    dualBoundaryLength (finRotate n) (fun i => ℓ i) j = 1 / 2 * (ℓ ((j + n - 1) % n) + ℓ j) := by
  simp only [dualBoundaryLength, finRotate_symm_val]

end Dual

/-! ## Exterior angles (`close_exterior_angles`, `flatten_to_disk`) -/

section Angles

variable {β : Type*} [Fintype β]

/-- The defect `2π − Σ k` of a list of exterior angles: what `close_exterior_angles` distributes
and returns. (`ClosureReport::turning_error` is its negative.) -/
def exteriorAngleDefect (k : β → ℝ) : ℝ := 2 * π - ∑ j, k j

/-- The angles written by `close_exterior_angles`: `k_j + (d_j / Σ d)·(2π − Σ k)`, each angle
gets a share of the defect proportional to the dual length `d_j`. Rust:
`*k += l / total_length * defect`. -/
def closedExteriorAngle (k d : β → ℝ) (j : β) : ℝ :=
  k j + d j / (∑ i, d i) * exteriorAngleDefect k

/-- **After `close_exterior_angles` the exterior angles sum to `2π`**, for any weights `d` whose
sum is not zero (the code requires the sum to be positive).

[formula: bff-close-exterior-angles] -/
theorem sum_closedExteriorAngle (k : β → ℝ) {d : β → ℝ} (hd : ∑ i, d i ≠ 0) :
    ∑ j, closedExteriorAngle k d j = 2 * π := by
  simp only [closedExteriorAngle, Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.sum_div]
  rw [div_self hd, exteriorAngleDefect]
  ring

/-- If the returned defect is zero, `close_exterior_angles` changes nothing (whatever `d` is).

[formula: bff-close-exterior-angles] -/
theorem closedExteriorAngle_of_defect_eq_zero (d : β → ℝ) {k : β → ℝ}
    (hk : exteriorAngleDefect k = 0) : closedExteriorAngle k d = k := by
  funext j
  simp [closedExteriorAngle, hk]

/-- The defect of the closed angles is zero: applying `close_exterior_angles` twice is the same
as applying it once. -/
theorem exteriorAngleDefect_closedExteriorAngle (k : β → ℝ) {d : β → ℝ} (hd : ∑ i, d i ≠ 0) :
    exteriorAngleDefect (closedExteriorAngle k d) = 0 := by
  rw [exteriorAngleDefect, sum_closedExteriorAngle k hd, sub_self]

/-- The exterior angles of `flatten_to_disk`, proportional to the dual lengths:
`k_j = 2π d_j / Σ d`. Rust: `*k = 2.0 * PI * l / total`. -/
def uniformizationAngle (d : β → ℝ) (j : β) : ℝ := 2 * π * d j / ∑ i, d i

/-- **The exterior angles of `flatten_to_disk` sum to `2π`** if the dual lengths do not sum to
zero.

[formula: bff-uniformization-angles] -/
theorem sum_uniformizationAngle {d : β → ℝ} (hd : ∑ i, d i ≠ 0) :
    ∑ j, uniformizationAngle d j = 2 * π := by
  simp only [uniformizationAngle]
  rw [← Finset.sum_div, ← Finset.mul_sum, mul_div_assoc, div_self hd, mul_one]

/-- The exterior angles of `flatten_to_disk` are positive if the dual lengths are: the target
curve is strictly convex.

[formula: bff-uniformization-angles] -/
theorem uniformizationAngle_pos {d : β → ℝ} (hd : ∀ j, 0 < d j) (j : β) :
    0 < uniformizationAngle d j := by
  have hsum : 0 < ∑ i, d i :=
    Finset.sum_pos (fun i _ => hd i) ⟨j, Finset.mem_univ j⟩
  have := hd j
  unfold uniformizationAngle
  positivity

/-- The angles of `flatten_to_disk` are what `close_exterior_angles` makes of the zero angles. -/
theorem uniformizationAngle_eq_closedExteriorAngle (d : β → ℝ) (j : β) :
    uniformizationAngle d j = closedExteriorAngle 0 d j := by
  simp [uniformizationAngle, closedExteriorAngle, exteriorAngleDefect]
  ring

/-- With `d` the dual lengths of edge lengths `ℓ`: the angles sum to `2π` if `Σ ℓ ≠ 0`. -/
theorem sum_uniformizationAngle_dual (next : Equiv.Perm β) {ℓ : β → ℝ} (hℓ : ∑ i, ℓ i ≠ 0) :
    ∑ j, uniformizationAngle (dualBoundaryLength next ℓ) j = 2 * π :=
  sum_uniformizationAngle (by rwa [sum_dualBoundaryLength])

/-- With `d` the dual lengths of positive edge lengths `ℓ`: the angles are positive. -/
theorem uniformizationAngle_dual_pos (next : Equiv.Perm β) {ℓ : β → ℝ} (hℓ : ∀ i, 0 < ℓ i)
    (j : β) : 0 < uniformizationAngle (dualBoundaryLength next ℓ) j :=
  uniformizationAngle_pos (dualBoundaryLength_pos next hℓ) j

/-- **One iteration of `flatten_to_disk`**, as coded: from scale factors `u` and positive input
lengths `ℓ`, the angles `2π d_j / Σ d` with `d` the dual lengths of the target lengths
`e^{(u_i + u_{next i})/2} ℓ_i` sum to `2π` and are all positive.

[formula: bff-uniformization-angles] -/
theorem flattenToDisk_angles [Nonempty β] (next : Equiv.Perm β) (u : β → ℝ) {ℓ : β → ℝ}
    (hℓ : ∀ i, 0 < ℓ i) :
    (∑ j, uniformizationAngle (dualBoundaryLength next (targetBoundaryLength next u ℓ)) j
        = 2 * π) ∧
      ∀ j, 0 < uniformizationAngle (dualBoundaryLength next (targetBoundaryLength next u ℓ)) j := by
  have hpos : ∀ i, 0 < targetBoundaryLength next u ℓ i :=
    fun i => targetBoundaryLength_pos next u (hℓ i)
  refine ⟨sum_uniformizationAngle_dual next ?_, uniformizationAngle_dual_pos next hpos⟩
  exact (Finset.sum_pos (fun i _ => hpos i) Finset.univ_nonempty).ne'

end Angles

/-! ## Tangents (`boundary_tangents`) -/

section Tangents

/-- The direction angle of edge `i`: `φ_0 = 0`, `φ_i = φ_{i−1} + k_i`. Rust: `phi` starts at `0`
and `phi += k` for `i > 0`; the angle `k_0` is not read. -/
def tangentAngle (k : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | i + 1 => tangentAngle k i + k (i + 1)

@[simp] theorem tangentAngle_zero (k : ℕ → ℝ) : tangentAngle k 0 = 0 := rfl

theorem tangentAngle_succ (k : ℕ → ℝ) (i : ℕ) :
    tangentAngle k (i + 1) = tangentAngle k i + k (i + 1) := rfl

/-- `k_0 + φ_i = k_0 + k_1 + … + k_i`. -/
theorem tangentAngle_add (k : ℕ → ℝ) (i : ℕ) :
    k 0 + tangentAngle k i = ∑ j ∈ range (i + 1), k j := by
  induction i with
  | zero => simp
  | succ i ih => rw [tangentAngle_succ, ← add_assoc, ih, Finset.sum_range_succ _ (i + 1)]

/-- The unit tangent of edge `i`, `T_i = (cos φ_i, sin φ_i)`, as the complex number
`exp (i φ_i)`. -/
def boundaryTangent (k : ℕ → ℝ) (i : ℕ) : ℂ := Complex.exp (tangentAngle k i * I)

/-- The first coordinate of the tangent is `cos φ_i`, as in the Rust code. -/
theorem boundaryTangent_re (k : ℕ → ℝ) (i : ℕ) :
    (boundaryTangent k i).re = Real.cos (tangentAngle k i) := by
  simp [boundaryTangent, Complex.exp_ofReal_mul_I_re]

/-- The second coordinate of the tangent is `sin φ_i`, as in the Rust code. -/
theorem boundaryTangent_im (k : ℕ → ℝ) (i : ℕ) :
    (boundaryTangent k i).im = Real.sin (tangentAngle k i) := by
  simp [boundaryTangent, Complex.exp_ofReal_mul_I_im]

/-- The first edge points along `+x`. -/
theorem boundaryTangent_zero (k : ℕ → ℝ) : boundaryTangent k 0 = 1 := by
  simp [boundaryTangent]

/-- **The tangents have unit length.**

[formula: bff-boundary-tangents] -/
theorem norm_boundaryTangent (k : ℕ → ℝ) (i : ℕ) : ‖boundaryTangent k i‖ = 1 :=
  Complex.norm_exp_ofReal_mul_I _

/-- **The turn at vertex `i + 1`, from edge `i` to edge `i + 1`, is the exterior angle
`k_{i+1}`**: `T_{i+1} = exp (i k_{i+1}) · T_i`.

[formula: bff-boundary-tangents] -/
theorem boundaryTangent_succ (k : ℕ → ℝ) (i : ℕ) :
    boundaryTangent k (i + 1) = Complex.exp (k (i + 1) * I) * boundaryTangent k i := by
  simp only [boundaryTangent, tangentAngle_succ]
  rw [← Complex.exp_add]
  push_cast
  ring_nf

/-- **The tangents close up exactly when the exterior angles sum to a multiple of `2π`**: for a
loop with `n` vertices, the turn at vertex `0`, from the last edge `n − 1` to the first edge, is
`k_0` if and only if `Σ_{j<n} k_j ∈ 2πℤ`.

[formula: bff-boundary-tangents] -/
theorem boundaryTangent_closes_iff (k : ℕ → ℝ) {n : ℕ} (hn : 0 < n) :
    boundaryTangent k 0 = Complex.exp (k 0 * I) * boundaryTangent k (n - 1) ↔
      ∃ m : ℤ, ∑ j ∈ range n, k j = m * (2 * π) := by
  obtain ⟨n, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hn.ne'
  rw [Nat.succ_sub_one, boundaryTangent_zero, boundaryTangent, ← Complex.exp_add, ← add_mul,
    ← Complex.ofReal_add, tangentAngle_add, eq_comm, Complex.exp_eq_one_iff]
  refine exists_congr fun m => ?_
  constructor
  · intro h
    have h' := congrArg Complex.im h
    simpa using h'
  · intro h
    rw [h]
    push_cast
    ring

/-- **If the exterior angles sum to `2π`, the turn at vertex `0` is `k_0` too.**

[formula: bff-boundary-tangents] -/
theorem boundaryTangent_closes (k : ℕ → ℝ) {n : ℕ} (hk : ∑ j ∈ range n, k j = 2 * π) :
    boundaryTangent k 0 = Complex.exp (k 0 * I) * boundaryTangent k (n - 1) := by
  have hn : 0 < n := by
    by_contra h0
    have hn0 : n = 0 := by omega
    rw [hn0, Finset.range_zero, Finset.sum_empty] at hk
    have := Real.pi_pos
    linarith
  exact (boundaryTangent_closes_iff k hn).2 ⟨1, by simpa using hk⟩

/-- The same for all vertices at once, in the cyclic indexing of the Rust code: if
`Σ_{j<n} k_j ∈ 2πℤ`, then at every vertex `i < n` the turn from edge `(i + n − 1) % n` to edge `i`
is `k_i`. -/
theorem boundaryTangent_turn (k : ℕ → ℝ) {n : ℕ}
    (hk : ∃ m : ℤ, ∑ j ∈ range n, k j = m * (2 * π)) {i : ℕ} (hi : i < n) :
    boundaryTangent k i = Complex.exp (k i * I) * boundaryTangent k ((i + n - 1) % n) := by
  rcases Nat.eq_zero_or_pos i with rfl | hpos
  · rw [Nat.zero_add, Nat.mod_eq_of_lt (by omega)]
    exact (boundaryTangent_closes_iff k (by omega)).2 hk
  · obtain ⟨i, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hpos.ne'
    have e : i + 1 + n - 1 = i + n := by omega
    rw [e, Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
    exact boundaryTangent_succ k i

/-- The tangents of the square (`k_j = π/2`) are `1, i, −1, −i, …`. -/
theorem boundaryTangent_square (i : ℕ) : boundaryTangent (fun _ => π / 2) i = I ^ i := by
  induction i with
  | zero => simp [boundaryTangent_zero]
  | succ i ih =>
    rw [boundaryTangent_succ, ih, pow_succ']
    congr 1
    have h : ((π / 2 : ℝ) : ℂ) = (π : ℂ) / 2 := by push_cast; rfl
    rw [h, Complex.exp_pi_div_two_mul_I]

end Tangents

/-! ## The curve (`integrate_curve`) and the spreading of the gap (`best_fit_curves`) -/

section Curve

/-- Vertex `p` of the boundary polygon with edge lengths `ℓ` and tangents `T`:
`γ_p = Σ_{i<p} ℓ_i T_i`, so `γ_0 = 0`. `integrate_curve` returns `γ_0, …, γ_{n−1}` and the end
point `γ_n`. -/
def boundaryCurve (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (p : ℕ) : ℂ := ∑ i ∈ range p, (ℓ i : ℂ) * T i

/-- The curve starts at the origin. -/
@[simp] theorem boundaryCurve_zero (ℓ : ℕ → ℝ) (T : ℕ → ℂ) : boundaryCurve ℓ T 0 = 0 := by
  simp [boundaryCurve]

/-- The loop of `integrate_curve`: `position += l * t`. -/
theorem boundaryCurve_succ (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (p : ℕ) :
    boundaryCurve ℓ T (p + 1) = boundaryCurve ℓ T p + ℓ p * T p :=
  Finset.sum_range_succ _ p

/-- First coordinate, as accumulated by the Rust code. -/
theorem boundaryCurve_re (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (p : ℕ) :
    (boundaryCurve ℓ T p).re = ∑ i ∈ range p, ℓ i * (T i).re := by
  simp [boundaryCurve]

/-- Second coordinate, as accumulated by the Rust code. -/
theorem boundaryCurve_im (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (p : ℕ) :
    (boundaryCurve ℓ T p).im = ∑ i ∈ range p, ℓ i * (T i).im := by
  simp [boundaryCurve]

/-- **Edge `p` of the polygon is `ℓ_p T_p`.**

[formula: bff-boundary-curve] -/
theorem boundaryCurve_succ_sub (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (p : ℕ) :
    boundaryCurve ℓ T (p + 1) - boundaryCurve ℓ T p = ℓ p * T p := by
  rw [boundaryCurve_succ, add_sub_cancel_left]

/-- **Edge `p` of the polygon has length `|ℓ_p|`** if `T_p` is a unit vector.

[formula: bff-boundary-curve] -/
theorem norm_boundaryCurve_succ_sub (ℓ : ℕ → ℝ) {T : ℕ → ℂ} {p : ℕ} (hT : ‖T p‖ = 1) :
    ‖boundaryCurve ℓ T (p + 1) - boundaryCurve ℓ T p‖ = |ℓ p| := by
  rw [boundaryCurve_succ_sub, norm_mul, hT, mul_one, Complex.norm_real, Real.norm_eq_abs]

/-- With the tangents of `boundary_tangents`, edge `p` has length `|ℓ_p|`. -/
theorem norm_boundaryCurve_boundaryTangent_succ_sub (ℓ k : ℕ → ℝ) (p : ℕ) :
    ‖boundaryCurve ℓ (boundaryTangent k) (p + 1) - boundaryCurve ℓ (boundaryTangent k) p‖
      = |ℓ p| :=
  norm_boundaryCurve_succ_sub ℓ (norm_boundaryTangent k p)

/-- The end point returned by `integrate_curve` is `γ_n = Σ_{i<n} ℓ_i T_i`. -/
theorem boundaryCurve_end (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (n : ℕ) :
    boundaryCurve ℓ T n = ∑ i ∈ range n, (ℓ i : ℂ) * T i := rfl

/-- **The polygon is closed iff `Σ_{i<n} ℓ_i T_i = 0`**: the closing edge, from the last vertex
`γ_{n−1}` back to `γ_0`, is the prescribed `ℓ_{n−1} T_{n−1}` if and only if the end point `γ_n`
is the origin.

[formula: bff-boundary-curve] -/
theorem boundaryCurve_closingEdge_iff (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (hn : 0 < n) :
    boundaryCurve ℓ T 0 - boundaryCurve ℓ T (n - 1) = ℓ (n - 1) * T (n - 1) ↔
      ∑ i ∈ range n, (ℓ i : ℂ) * T i = 0 := by
  obtain ⟨n, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hn.ne'
  rw [Nat.succ_sub_one, ← boundaryCurve_end, boundaryCurve_succ, boundaryCurve_zero]
  constructor
  · intro h; rw [← h]; ring
  · intro h; linear_combination -h

/-- Arc length up to vertex `j`: `arc_j = Σ_{i<j} ℓ_i`; `boundaryArc ℓ n` is the total length. -/
def boundaryArc (ℓ : ℕ → ℝ) (j : ℕ) : ℝ := ∑ i ∈ range j, ℓ i

/-- The polygon after the gap `γ_n` is spread linearly along the loop:
`p'_j = γ_j − (arc_j / total)·γ_n`. The code returns `p'_0, …, p'_{n−1}`. -/
def spreadCurve (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (n j : ℕ) : ℂ :=
  boundaryCurve ℓ T j - ((boundaryArc ℓ j / boundaryArc ℓ n : ℝ) : ℂ) * boundaryCurve ℓ T n

/-- The same with the guard of the Rust code:
`let s = if total > 0.0 { arc / total } else { 0.0 }; p -= s * end`. -/
def spreadCurveCode (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (n j : ℕ) : ℂ :=
  boundaryCurve ℓ T j -
    ((if 0 < boundaryArc ℓ n then boundaryArc ℓ j / boundaryArc ℓ n else 0 : ℝ) : ℂ) *
      boundaryCurve ℓ T n

/-- For a positive total length the code evaluates `spreadCurve`. -/
theorem spreadCurveCode_of_pos (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (h : 0 < boundaryArc ℓ n)
    (j : ℕ) : spreadCurveCode ℓ T n j = spreadCurve ℓ T n j := by
  simp [spreadCurveCode, spreadCurve, h]

/-- If the total length is not positive, the code leaves the polygon as it is: the gap is then
not removed. -/
theorem spreadCurveCode_of_not_pos (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (h : ¬ 0 < boundaryArc ℓ n)
    (j : ℕ) : spreadCurveCode ℓ T n j = boundaryCurve ℓ T j := by
  simp [spreadCurveCode, h]

/-- The corrected polygon starts at the origin. -/
theorem spreadCurve_zero (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (n : ℕ) : spreadCurve ℓ T n 0 = 0 := by
  simp [spreadCurve, boundaryArc]

/-- **The corrected polygon is exactly closed**: the formula at `j = n` gives the origin again,
if the total length is not zero.

[formula: bff-boundary-curve] -/
theorem spreadCurve_end (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (h : boundaryArc ℓ n ≠ 0) :
    spreadCurve ℓ T n n = 0 := by
  simp [spreadCurve, div_self h]

/-- **Each edge changes by its share of the gap**: edge `j` of the corrected polygon is
`ℓ_j T_j − (ℓ_j / total)·γ_n`. (No hypothesis on the total length is needed for this identity.)

[formula: bff-boundary-curve] -/
theorem spreadCurve_succ_sub (ℓ : ℕ → ℝ) (T : ℕ → ℂ) (n j : ℕ) :
    spreadCurve ℓ T n (j + 1) - spreadCurve ℓ T n j =
      ℓ j * T j - ((ℓ j / boundaryArc ℓ n : ℝ) : ℂ) * boundaryCurve ℓ T n := by
  have harc : boundaryArc ℓ (j + 1) = boundaryArc ℓ j + ℓ j := Finset.sum_range_succ _ j
  simp only [spreadCurve, boundaryCurve_succ ℓ T j, harc]
  push_cast
  ring

/-- The closing edge of the corrected polygon, from its last vertex `p'_{n−1}` back to `p'_0`,
follows the same formula as the other edges. -/
theorem spreadCurve_closingEdge (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (hn : 0 < n)
    (h : boundaryArc ℓ n ≠ 0) :
    spreadCurve ℓ T n 0 - spreadCurve ℓ T n (n - 1) =
      ℓ (n - 1) * T (n - 1) - ((ℓ (n - 1) / boundaryArc ℓ n : ℝ) : ℂ) * boundaryCurve ℓ T n := by
  have h1 := spreadCurve_succ_sub ℓ T n (n - 1)
  rw [Nat.sub_add_cancel hn, spreadCurve_end ℓ T h] at h1
  rw [spreadCurve_zero, ← h1]

/-- A polygon that is already closed is not changed. -/
theorem spreadCurve_of_closed (ℓ : ℕ → ℝ) (T : ℕ → ℂ) {n : ℕ} (h : boundaryCurve ℓ T n = 0)
    (j : ℕ) : spreadCurve ℓ T n j = boundaryCurve ℓ T j := by
  simp [spreadCurve, h]

/-! ### The unit square

`n = 4`, `k_j = π/2`, `ℓ_j = 1` (the test `a_closed_curve_is_reproduced` of `curve.rs`). -/

/-- The exterior angles of the square sum to `2π`. -/
example : ∑ j ∈ range 4, (fun _ : ℕ => π / 2) j = 2 * π := by
  simp; ring

/-- The tangents of the square close up. -/
example :
    boundaryTangent (fun _ => π / 2) 0 =
      Complex.exp ((π / 2 : ℝ) * I) * boundaryTangent (fun _ => π / 2) (4 - 1) :=
  boundaryTangent_closes _ (by simp; ring)

/-- The vertices of the square are `0, 1, 1 + i, i`. -/
example :
    boundaryCurve (fun _ => 1) (boundaryTangent fun _ => π / 2) 1 = 1 ∧
    boundaryCurve (fun _ => 1) (boundaryTangent fun _ => π / 2) 2 = 1 + I ∧
    boundaryCurve (fun _ => 1) (boundaryTangent fun _ => π / 2) 3 = I := by
  simp [boundaryCurve, Finset.sum_range_succ, boundaryTangent_square, pow_succ]

/-- The square closes: its end point `γ_4` is the origin. -/
example : boundaryCurve (fun _ => 1) (boundaryTangent fun _ => π / 2) 4 = 0 := by
  simp [boundaryCurve, Finset.sum_range_succ, boundaryTangent_square, pow_succ]

end Curve

/-! ## Closing by lengths (`best_fit_curves`, `close_group`; paper Eq. 19–20)

The closing problem is

```text
minimize  ½ Σ_d N_d (ℓ̃_d − ℓ*_d)²   subject to   T ℓ̃ = 0,
```

with `δ` the length degrees of freedom, `κ` the closing constraints (in the code two per loop: the
row `2·block + axis`), `T : Matrix κ δ ℝ` the constraint matrix, `N_d > 0` the weights and `ℓ*`
the target lengths.

Dictionary with `close_group`:

* `dofs` enumerates `δ`; `column(dof)` returns the non-zero entries `(rows, values)` of column `d`
  of `T`: the two components of the tangent of each side of the degree of freedom, added up when
  both sides lie in the same loop;
* `dof.weight` is `1 / N_d` (`N_d = Σ 1/ℓ_i` over the sides), `dof.target` is `ℓ*_d`;
* `s[rows[a] * m + rows[b]] += dof.weight * values[a] * values[b]` assembles
  `S = closingMatrix T N = T · diag(1/N) · Tᵀ` (`closingMatrix_eq`);
* `r[rows[a]] += dof.target * values[a]` assembles `r = T *ᵥ ℓ*`;
* `lambda` is the solution `lam` of `S lam = r` returned by `solve_semidefinite`;
* `new = dof.target - dof.weight * dot`, with `dot = Σ_a values[a] * lambda[rows[a]]`, is
  `closedLengths T N ℓ* lam d = ℓ*_d − (1/N_d)·(Tᵀ lam)_d` (`closedLengths_eq`).

The factor `½` of the objective does not change the minimizer and is left out below. -/

section Closing

variable {κ δ : Type*} [Fintype κ] [Fintype δ]

/-- The matrix `S = T · diag(1/N) · Tᵀ` of the closing system, entry by entry as the code
accumulates it: `S_ab = Σ_d (1/N_d)·T_ad·T_bd`. -/
def closingMatrix (T : Matrix κ δ ℝ) (N : δ → ℝ) : Matrix κ κ ℝ :=
  Matrix.of fun a b => ∑ d, 1 / N d * T a d * T b d

/-- The adjusted lengths, as the code computes them from the multipliers `lam`:
`ℓ̃_d = ℓ*_d − (1/N_d)·Σ_a T_ad·lam_a`. -/
def closedLengths (T : Matrix κ δ ℝ) (N ℓstar : δ → ℝ) (lam : κ → ℝ) (d : δ) : ℝ :=
  ℓstar d - 1 / N d * ∑ a, T a d * lam a

omit [Fintype κ] in
/-- The code's `S` is `T N⁻¹ Tᵀ`. -/
theorem closingMatrix_eq [DecidableEq δ] (T : Matrix κ δ ℝ) (N : δ → ℝ) :
    closingMatrix T N = T * Matrix.diagonal (fun d => 1 / N d) * Tᵀ := by
  ext a b
  rw [Matrix.mul_apply]
  simp only [closingMatrix, Matrix.of_apply, Matrix.mul_diagonal, Matrix.transpose_apply]
  exact Finset.sum_congr rfl fun d _ => by ring

/-- The code's update is `ℓ̃ = ℓ* − N⁻¹ Tᵀ lam`. -/
theorem closedLengths_eq [DecidableEq δ] (T : Matrix κ δ ℝ) (N ℓstar : δ → ℝ) (lam : κ → ℝ) :
    closedLengths T N ℓstar lam =
      ℓstar - Matrix.diagonal (fun d => 1 / N d) *ᵥ (Tᵀ *ᵥ lam) := by
  funext d
  rw [closedLengths, Pi.sub_apply, Matrix.mulVec_diagonal]
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]

/-- The constraint residual of the adjusted lengths: `T ℓ̃ = T ℓ* − S lam`. -/
theorem mulVec_closedLengths (T : Matrix κ δ ℝ) (N ℓstar : δ → ℝ) (lam : κ → ℝ) :
    T *ᵥ closedLengths T N ℓstar lam = T *ᵥ ℓstar - closingMatrix T N *ᵥ lam := by
  funext a
  simp only [Matrix.mulVec, dotProduct, closedLengths, closingMatrix, Matrix.of_apply,
    Pi.sub_apply, mul_sub, Finset.sum_sub_distrib, Finset.sum_mul, Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun b _ => by ring

/-- **The adjusted lengths satisfy the closing constraints exactly when the multipliers solve
`S lam = T ℓ*`.** No assumption on the weights.

[formula: bff-closing-lengths] -/
theorem closedLengths_closes_iff (T : Matrix κ δ ℝ) (N ℓstar : δ → ℝ) (lam : κ → ℝ) :
    T *ᵥ closedLengths T N ℓstar lam = 0 ↔ closingMatrix T N *ᵥ lam = T *ᵥ ℓstar := by
  rw [mulVec_closedLengths, sub_eq_zero, eq_comm]

/-- If `S lam = T ℓ*`, the adjusted lengths close: `T ℓ̃ = 0`. -/
theorem closedLengths_closes {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ} {lam : κ → ℝ}
    (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar) :
    T *ᵥ closedLengths T N ℓstar lam = 0 :=
  (closedLengths_closes_iff T N ℓstar lam).2 hlam

/-- Pythagoras for the closing problem: for lengths `ℓ'` that close, the cost of `ℓ'` is the cost
of the adjusted lengths `ℓ̃` plus the weighted distance from `ℓ'` to `ℓ̃`. -/
theorem closingCost_eq_add {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ} {lam : κ → ℝ}
    (hN : ∀ d, N d ≠ 0) (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar)
    {ℓ' : δ → ℝ} (hℓ' : T *ᵥ ℓ' = 0) :
    ∑ d, N d * (ℓ' d - ℓstar d) ^ 2 =
      ∑ d, N d * (closedLengths T N ℓstar lam d - ℓstar d) ^ 2 +
        ∑ d, N d * (ℓ' d - closedLengths T N ℓstar lam d) ^ 2 := by
  have hc := closedLengths_closes hlam
  set ℓt := closedLengths T N ℓstar lam with hℓt
  have hcross : ∑ d, (∑ a, T a d * lam a) * (ℓ' d - ℓt d) = 0 := by
    have e : ∑ d, (∑ a, T a d * lam a) * (ℓ' d - ℓt d) =
        ∑ a, lam a * ((T *ᵥ ℓ') a - (T *ᵥ ℓt) a) := by
      simp only [Matrix.mulVec, dotProduct, Finset.sum_mul, Finset.mul_sum,
        ← Finset.sum_sub_distrib]
      rw [Finset.sum_comm]
      exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun d _ => by ring
    rw [e, hℓ', hc]
    simp
  have hpt : ∀ d, N d * (ℓ' d - ℓstar d) ^ 2 =
      N d * (ℓt d - ℓstar d) ^ 2 + N d * (ℓ' d - ℓt d) ^ 2 -
        2 * ((∑ a, T a d * lam a) * (ℓ' d - ℓt d)) := by
    intro d
    have h : N d * (ℓt d - ℓstar d) = -(∑ a, T a d * lam a) := by
      rw [hℓt, closedLengths]
      field_simp [hN d]
      ring
    linear_combination (2 * (ℓ' d - ℓt d)) * h
  rw [Finset.sum_congr rfl fun d _ => hpt d, Finset.sum_sub_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, hcross]
  ring

/-- **The adjusted lengths are the closest lengths that close** (paper Eq. 19): with positive
weights `N` and multipliers solving `S lam = T ℓ*`, every `ℓ'` with `T ℓ' = 0` is at least as far
from the target lengths, in the norm `Σ_d N_d (·)²`, as `ℓ̃ = ℓ* − N⁻¹ Tᵀ lam`.

[formula: bff-closing-lengths] -/
theorem closedLengths_optimal {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ} {lam : κ → ℝ}
    (hN : ∀ d, 0 < N d) (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar)
    {ℓ' : δ → ℝ} (hℓ' : T *ᵥ ℓ' = 0) :
    ∑ d, N d * (closedLengths T N ℓstar lam d - ℓstar d) ^ 2 ≤
      ∑ d, N d * (ℓ' d - ℓstar d) ^ 2 := by
  rw [closingCost_eq_add (fun d => (hN d).ne') hlam hℓ']
  have : 0 ≤ ∑ d, N d * (ℓ' d - closedLengths T N ℓstar lam d) ^ 2 :=
    Finset.sum_nonneg fun d _ => mul_nonneg (hN d).le (sq_nonneg _)
  linarith

/-- **The minimizer is unique**: lengths `ℓ'` that close and are not farther from the target
lengths than the adjusted lengths (in particular: at the same distance) are the adjusted lengths.

[formula: bff-closing-lengths] -/
theorem closedLengths_unique {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ} {lam : κ → ℝ}
    (hN : ∀ d, 0 < N d) (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar)
    {ℓ' : δ → ℝ} (hℓ' : T *ᵥ ℓ' = 0)
    (hle : ∑ d, N d * (ℓ' d - ℓstar d) ^ 2 ≤
      ∑ d, N d * (closedLengths T N ℓstar lam d - ℓstar d) ^ 2) :
    ℓ' = closedLengths T N ℓstar lam := by
  rw [closingCost_eq_add (fun d => (hN d).ne') hlam hℓ'] at hle
  have hnn : ∀ d ∈ Finset.univ, 0 ≤ N d * (ℓ' d - closedLengths T N ℓstar lam d) ^ 2 :=
    fun d _ => mul_nonneg (hN d).le (sq_nonneg _)
  have hzero : ∑ d, N d * (ℓ' d - closedLengths T N ℓstar lam d) ^ 2 = 0 :=
    le_antisymm (by linarith) (Finset.sum_nonneg hnn)
  funext d
  have hd := (Finset.sum_eq_zero_iff_of_nonneg hnn).1 hzero d (Finset.mem_univ d)
  have hsq : (ℓ' d - closedLengths T N ℓstar lam d) ^ 2 = 0 := by
    rcases mul_eq_zero.1 hd with h | h
    · exact absurd h (hN d).ne'
    · exact h
  exact sub_eq_zero.1 ((pow_eq_zero_iff two_ne_zero).1 hsq)

/-- The adjusted lengths do not depend on which solution of `S lam = T ℓ*` is used (relevant
when `S` is singular). -/
theorem closedLengths_eq_of_solutions {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ} {lam lam' : κ → ℝ}
    (hN : ∀ d, 0 < N d) (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar)
    (hlam' : closingMatrix T N *ᵥ lam' = T *ᵥ ℓstar) :
    closedLengths T N ℓstar lam' = closedLengths T N ℓstar lam :=
  closedLengths_unique hN hlam (closedLengths_closes hlam')
    (closedLengths_optimal hN hlam' (closedLengths_closes hlam))

/-- For positive weights the closing system `S lam = T ℓ*` has a solution, whatever the rank of
`T` is: the right-hand side lies in the range of `S = T N⁻¹ Tᵀ`, which is the range of `T`. -/
theorem exists_closing_multiplier {T : Matrix κ δ ℝ} {N : δ → ℝ} (hN : ∀ d, 0 < N d)
    (ℓstar : δ → ℝ) : ∃ lam : κ → ℝ, closingMatrix T N *ᵥ lam = T *ᵥ ℓstar := by
  classical
  set R : Matrix δ δ ℝ := Matrix.diagonal fun d => (Real.sqrt (N d))⁻¹ with hR
  have hRR : R * Rᵀ = Matrix.diagonal fun d => 1 / N d := by
    rw [hR, Matrix.diagonal_transpose, Matrix.diagonal_mul_diagonal]
    congr 1
    funext d
    rw [← mul_inv, Real.mul_self_sqrt (hN d).le, one_div]
  have hS : closingMatrix T N = (T * R) * (T * R)ᵀ := by
    rw [closingMatrix_eq, Matrix.transpose_mul, ← hRR]
    simp only [Matrix.mul_assoc]
  have hdet : IsUnit R.det := by
    rw [hR, Matrix.det_diagonal, isUnit_iff_ne_zero]
    exact Finset.prod_ne_zero_iff.2 fun d _ => inv_ne_zero (Real.sqrt_ne_zero'.2 (hN d))
  have hrank : (closingMatrix T N).rank = T.rank := by
    rw [hS, Matrix.rank_self_mul_transpose, Matrix.rank_mul_eq_left_of_isUnit_det R T hdet]
  have hle : LinearMap.range (closingMatrix T N).mulVecLin ≤ LinearMap.range T.mulVecLin := by
    rw [closingMatrix_eq, Matrix.mul_assoc, Matrix.mulVecLin_mul]
    exact LinearMap.range_comp_le_range _ _
  have heq := Submodule.eq_of_le_of_finrank_eq hle hrank
  have hmem : T *ᵥ ℓstar ∈ LinearMap.range T.mulVecLin := ⟨ℓstar, rfl⟩
  rw [← heq] at hmem
  obtain ⟨lam, hlam⟩ := hmem
  exact ⟨lam, hlam⟩

/-- If `S` is invertible, the multipliers are `lam = S⁻¹ T ℓ*`. -/
theorem closing_multiplier_eq [DecidableEq κ] {T : Matrix κ δ ℝ} {N ℓstar : δ → ℝ}
    {lam : κ → ℝ} (hS : IsUnit (closingMatrix T N).det)
    (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar) :
    lam = (closingMatrix T N)⁻¹ *ᵥ (T *ᵥ ℓstar) := by
  rw [← hlam, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ hS, Matrix.one_mulVec]

/-- If `S` is invertible, `lam = S⁻¹ T ℓ*` solves the closing system. -/
theorem closing_multiplier_solves [DecidableEq κ] {T : Matrix κ δ ℝ} {N : δ → ℝ}
    (hS : IsUnit (closingMatrix T N).det) (ℓstar : δ → ℝ) :
    closingMatrix T N *ᵥ ((closingMatrix T N)⁻¹ *ᵥ (T *ᵥ ℓstar)) = T *ᵥ ℓstar := by
  rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hS, Matrix.one_mulVec]

/-- The formula of the paper (Eq. 20): `ℓ̃ = ℓ* − N⁻¹ Tᵀ (T N⁻¹ Tᵀ)⁻¹ T ℓ*`, with
`N⁻¹ = diag(1/N_d)`. -/
def closedLengthsPaper [DecidableEq κ] [DecidableEq δ] (T : Matrix κ δ ℝ) (N ℓstar : δ → ℝ) :
    δ → ℝ :=
  ℓstar - (Matrix.diagonal (fun d => 1 / N d) * Tᵀ *
    (T * Matrix.diagonal (fun d => 1 / N d) * Tᵀ)⁻¹ * T) *ᵥ ℓstar

/-- The formula of the paper is the code's update with `lam = S⁻¹ T ℓ*`. -/
theorem closedLengthsPaper_eq [DecidableEq κ] [DecidableEq δ] (T : Matrix κ δ ℝ)
    (N ℓstar : δ → ℝ) :
    closedLengthsPaper T N ℓstar =
      closedLengths T N ℓstar ((closingMatrix T N)⁻¹ *ᵥ (T *ᵥ ℓstar)) := by
  rw [closedLengths_eq, closedLengthsPaper, closingMatrix_eq]
  simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]

/-- If `S` is invertible, the code's lengths are those of the paper's formula. -/
theorem closedLengths_eq_paper [DecidableEq κ] [DecidableEq δ] {T : Matrix κ δ ℝ}
    {N ℓstar : δ → ℝ} {lam : κ → ℝ} (hS : IsUnit (closingMatrix T N).det)
    (hlam : closingMatrix T N *ᵥ lam = T *ᵥ ℓstar) :
    closedLengths T N ℓstar lam = closedLengthsPaper T N ℓstar := by
  rw [closedLengthsPaper_eq, ← closing_multiplier_eq hS hlam]

/-- **The formula of the paper (Eq. 20) solves the closing problem (Eq. 19)**: for positive
weights and invertible `S = T N⁻¹ Tᵀ`, the lengths `ℓ̃ = ℓ* − N⁻¹ Tᵀ (T N⁻¹ Tᵀ)⁻¹ T ℓ*` close,
no lengths that close are nearer to the target lengths, and lengths that close and are as near
are equal to `ℓ̃`.

[formula: bff-closing-lengths] -/
theorem closedLengthsPaper_isMinimizer [DecidableEq κ] [DecidableEq δ] {T : Matrix κ δ ℝ}
    {N : δ → ℝ} (hN : ∀ d, 0 < N d) (hS : IsUnit (closingMatrix T N).det) (ℓstar : δ → ℝ) :
    T *ᵥ closedLengthsPaper T N ℓstar = 0 ∧
      ∀ ℓ' : δ → ℝ, T *ᵥ ℓ' = 0 →
        ∑ d, N d * (closedLengthsPaper T N ℓstar d - ℓstar d) ^ 2 ≤
            ∑ d, N d * (ℓ' d - ℓstar d) ^ 2 ∧
          (∑ d, N d * (ℓ' d - ℓstar d) ^ 2 ≤
              ∑ d, N d * (closedLengthsPaper T N ℓstar d - ℓstar d) ^ 2 →
            ℓ' = closedLengthsPaper T N ℓstar) := by
  have hlam := closing_multiplier_solves hS ℓstar
  rw [closedLengthsPaper_eq]
  exact ⟨closedLengths_closes hlam, fun ℓ' hℓ' =>
    ⟨closedLengths_optimal hN hlam hℓ', closedLengths_unique hN hlam hℓ'⟩⟩

/-- A tiny instance of the hypotheses: one constraint `ℓ̃_0 = ℓ̃_1`, unit weights, targets `1`
and `3`. The multiplier `lam = −1` solves `S lam = T ℓ*` (`2·(−1) = 1 − 3`) and the adjusted
lengths are `2` and `2`. -/
example :
    closingMatrix !![(1 : ℝ), -1] ![1, 1] *ᵥ ![-1] = !![(1 : ℝ), -1] *ᵥ ![1, 3] ∧
      closedLengths !![(1 : ℝ), -1] ![1, 1] ![1, 3] ![-1] = ![2, 2] := by
  constructor
  · funext a
    fin_cases a
    simp [closingMatrix, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    norm_num
  · funext d
    fin_cases d <;> simp [closedLengths] <;> norm_num

/-! ### One loop without shared lengths

`κ = Fin 2` (the two coordinates), `δ = Fin n` (the edges), column `d` of `T` is the tangent
`T_d`, and `N_d = 1/ℓ_d` with `ℓ` the edge lengths of the input metric, so `1/N_d = ℓ_d`. The
sequences are indexed by `ℕ` as in the sections on tangents and curves; only the entries
`d < n` matter. -/

/-- The `2 × n` constraint matrix of one loop: row `0` holds the first coordinates of the
tangents, row `1` the second coordinates. -/
def loopMatrix (n : ℕ) (Tan : ℕ → ℂ) : Matrix (Fin 2) (Fin n) ℝ :=
  Matrix.of ![fun d => (Tan d).re, fun d => (Tan d).im]

/-- Row `0` of `T x` is the first coordinate of the end point `Σ_{d<n} x_d T_d`. -/
theorem loopMatrix_mulVec_zero (n : ℕ) (Tan : ℕ → ℂ) (x : ℕ → ℝ) :
    (loopMatrix n Tan *ᵥ fun d : Fin n => x d) 0 = (boundaryCurve x Tan n).re := by
  rw [boundaryCurve_re, ← Fin.sum_univ_eq_sum_range (fun d => x d * (Tan d).re) n]
  simp [loopMatrix, Matrix.mulVec, dotProduct, mul_comm]

/-- Row `1` of `T x` is the second coordinate of the end point `Σ_{d<n} x_d T_d`. -/
theorem loopMatrix_mulVec_one (n : ℕ) (Tan : ℕ → ℂ) (x : ℕ → ℝ) :
    (loopMatrix n Tan *ᵥ fun d : Fin n => x d) 1 = (boundaryCurve x Tan n).im := by
  rw [boundaryCurve_im, ← Fin.sum_univ_eq_sum_range (fun d => x d * (Tan d).im) n]
  simp [loopMatrix, Matrix.mulVec, dotProduct, mul_comm]

/-- The closing constraint `T x = 0` of one loop says that the polygon with edge lengths `x`
closes. -/
theorem loopMatrix_mulVec_eq_zero_iff (n : ℕ) (Tan : ℕ → ℂ) (x : ℕ → ℝ) :
    (loopMatrix n Tan *ᵥ fun d : Fin n => x d) = 0 ↔ boundaryCurve x Tan n = 0 := by
  rw [funext_iff, Fin.forall_fin_two, loopMatrix_mulVec_zero, loopMatrix_mulVec_one,
    Complex.ext_iff]
  rfl

/-- The adjusted lengths of one loop: `ℓ̃_d = ℓ*_d − ℓ_d·⟨T_d, lam⟩` (as written out in the test
`matches_the_formula_of_the_paper_for_one_loop` of `curve.rs`). -/
def loopClosedLength (ℓ ℓstar : ℕ → ℝ) (Tan : ℕ → ℂ) (lam : Fin 2 → ℝ) (d : ℕ) : ℝ :=
  ℓstar d - ℓ d * ((Tan d).re * lam 0 + (Tan d).im * lam 1)

/-- `loopClosedLength` is `closedLengths` for the matrix of one loop and `N_d = 1/ℓ_d`. -/
theorem loopClosedLength_eq (n : ℕ) (ℓ ℓstar : ℕ → ℝ) (Tan : ℕ → ℂ) (lam : Fin 2 → ℝ) :
    (fun d : Fin n => loopClosedLength ℓ ℓstar Tan lam d) =
      closedLengths (loopMatrix n Tan) (fun d => 1 / ℓ d) (fun d => ℓstar d) lam := by
  funext d
  simp [loopClosedLength, closedLengths, loopMatrix, Fin.sum_univ_two]

/-- The closing system `S lam = T ℓ*` of one loop in vector form:
`Σ_{d<n} ℓ_d ⟨T_d, lam⟩ T_d = Σ_{d<n} ℓ*_d T_d`, the right-hand side being the end point of the
curve with the target lengths. -/
theorem loopClosing_system_iff (n : ℕ) (ℓ ℓstar : ℕ → ℝ) (Tan : ℕ → ℂ) (lam : Fin 2 → ℝ) :
    closingMatrix (loopMatrix n Tan) (fun d => 1 / ℓ d) *ᵥ lam =
        loopMatrix n Tan *ᵥ (fun d => ℓstar d) ↔
      ∑ d ∈ range n, ((ℓ d * ((Tan d).re * lam 0 + (Tan d).im * lam 1) : ℝ) : ℂ) * Tan d =
        boundaryCurve ℓstar Tan n := by
  rw [← closedLengths_closes_iff, ← loopClosedLength_eq, loopMatrix_mulVec_eq_zero_iff]
  simp only [boundaryCurve, loopClosedLength]
  push_cast
  simp only [sub_mul, Finset.sum_sub_distrib, sub_eq_zero]
  exact eq_comm

/-- **The polygon with the adjusted lengths closes**: for one loop, if `lam` solves the closing
system, then `Σ_{d<n} ℓ̃_d T_d = 0`.

[formula: bff-closing-lengths] -/
theorem boundaryCurve_loopClosedLength (n : ℕ) {ℓ ℓstar : ℕ → ℝ} {Tan : ℕ → ℂ}
    {lam : Fin 2 → ℝ}
    (hlam : closingMatrix (loopMatrix n Tan) (fun d => 1 / ℓ d) *ᵥ lam =
      loopMatrix n Tan *ᵥ (fun d => ℓstar d)) :
    boundaryCurve (loopClosedLength ℓ ℓstar Tan lam) Tan n = 0 := by
  rw [← loopMatrix_mulVec_eq_zero_iff, loopClosedLength_eq]
  exact closedLengths_closes hlam

/-- After closing by lengths there is no gap left to spread: `spreadCurve` returns the polygon
of the adjusted lengths. -/
theorem spreadCurve_loopClosedLength (n : ℕ) {ℓ ℓstar : ℕ → ℝ} {Tan : ℕ → ℂ} {lam : Fin 2 → ℝ}
    (hlam : closingMatrix (loopMatrix n Tan) (fun d => 1 / ℓ d) *ᵥ lam =
      loopMatrix n Tan *ᵥ (fun d => ℓstar d)) (j : ℕ) :
    spreadCurve (loopClosedLength ℓ ℓstar Tan lam) Tan n j =
      boundaryCurve (loopClosedLength ℓ ℓstar Tan lam) Tan j :=
  spreadCurve_of_closed _ _ (boundaryCurve_loopClosedLength n hlam) j

/-- **The adjusted lengths of one loop are the closest lengths of a closed polygon with the
given tangents** (paper Eq. 19): for positive input lengths `ℓ_d`, every `ℓ'` whose polygon closes
has `Σ_{d<n} (ℓ̃_d − ℓ*_d)²/ℓ_d ≤ Σ_{d<n} (ℓ'_d − ℓ*_d)²/ℓ_d`, and if `ℓ'` is as near to the
target lengths then `ℓ'_d = ℓ̃_d` for all `d < n`.

[formula: bff-closing-lengths] -/
theorem loopClosedLength_optimal (n : ℕ) {ℓ ℓstar : ℕ → ℝ} {Tan : ℕ → ℂ} {lam : Fin 2 → ℝ}
    (hℓ : ∀ d < n, 0 < ℓ d)
    (hlam : closingMatrix (loopMatrix n Tan) (fun d => 1 / ℓ d) *ᵥ lam =
      loopMatrix n Tan *ᵥ (fun d => ℓstar d))
    {ℓ' : ℕ → ℝ} (hℓ' : boundaryCurve ℓ' Tan n = 0) :
    ∑ d ∈ range n, (loopClosedLength ℓ ℓstar Tan lam d - ℓstar d) ^ 2 / ℓ d ≤
        ∑ d ∈ range n, (ℓ' d - ℓstar d) ^ 2 / ℓ d ∧
      (∑ d ∈ range n, (ℓ' d - ℓstar d) ^ 2 / ℓ d ≤
          ∑ d ∈ range n, (loopClosedLength ℓ ℓstar Tan lam d - ℓstar d) ^ 2 / ℓ d →
        ∀ d < n, ℓ' d = loopClosedLength ℓ ℓstar Tan lam d) := by
  have hN : ∀ d : Fin n, 0 < 1 / ℓ d := fun d => one_div_pos.2 (hℓ d d.isLt)
  have hc : loopMatrix n Tan *ᵥ (fun d : Fin n => ℓ' d) = 0 :=
    (loopMatrix_mulVec_eq_zero_iff n Tan ℓ').2 hℓ'
  have e1 : ∑ d ∈ range n, (loopClosedLength ℓ ℓstar Tan lam d - ℓstar d) ^ 2 / ℓ d =
      ∑ d : Fin n, 1 / ℓ d * (closedLengths (loopMatrix n Tan) (fun d => 1 / ℓ d)
        (fun d => ℓstar d) lam d - ℓstar d) ^ 2 := by
    rw [← loopClosedLength_eq,
      ← Fin.sum_univ_eq_sum_range
        (fun d => (loopClosedLength ℓ ℓstar Tan lam d - ℓstar d) ^ 2 / ℓ d) n]
    exact Finset.sum_congr rfl fun d _ => by ring
  have e2 : ∑ d ∈ range n, (ℓ' d - ℓstar d) ^ 2 / ℓ d =
      ∑ d : Fin n, 1 / ℓ d * (ℓ' d - ℓstar d) ^ 2 := by
    rw [← Fin.sum_univ_eq_sum_range (fun d => (ℓ' d - ℓstar d) ^ 2 / ℓ d) n]
    exact Finset.sum_congr rfl fun d _ => by ring
  rw [e1, e2]
  refine ⟨closedLengths_optimal hN hlam hc, fun hle d hd => ?_⟩
  have h := closedLengths_unique hN hlam hc hle
  rw [← loopClosedLength_eq] at h
  exact congrFun h ⟨d, hd⟩

/-! ### A rectangle whose opposite sides disagree

The test `lengths_are_adjusted_minimally_to_close` of `curve.rs`: the tangents of the square,
input lengths `1`, target lengths `1, 2, 6/5, 2`. The multipliers `lam = (−1/10, 0)` solve the
closing system, and the adjusted lengths are `11/10, 2, 11/10, 2`: the two sides that disagree
are averaged. -/

/-- The hypotheses of the one-loop statements hold for this instance. -/
example :
    closingMatrix (loopMatrix 4 (boundaryTangent fun _ => π / 2)) (fun _ => 1 / 1) *ᵥ
        ![-1 / 10, 0] =
      loopMatrix 4 (boundaryTangent fun _ => π / 2) *ᵥ
        (fun d : Fin 4 => if (d : ℕ) = 0 then 1 else if (d : ℕ) = 2 then 6 / 5 else 2) := by
  rw [loopClosing_system_iff 4 (fun _ => 1)
    (fun d => if d = 0 then 1 else if d = 2 then 6 / 5 else 2)]
  simp [boundaryCurve, Finset.sum_range_succ, boundaryTangent_square, pow_succ]
  norm_num [Complex.ext_iff]

/-- The adjusted lengths of this instance. -/
example :
    (fun d : Fin 4 => loopClosedLength (fun _ => 1)
      (fun d => if d = 0 then 1 else if d = 2 then 6 / 5 else 2)
      (boundaryTangent fun _ => π / 2) ![-1 / 10, 0] d) = ![11 / 10, 2, 11 / 10, 2] := by
  funext d
  fin_cases d <;>
    simp [loopClosedLength, boundaryTangent_square, pow_succ] <;> norm_num

end Closing

/-! ## The right-hand side for the harmonic conjugate (`hilbert_transform`; paper Eq. 17)

The paper writes the conformal energy with the matrix `U` of its Eq. 16,
`aᵀUb = ½ Σ (a_j b_i − a_i b_j)` over the boundary edges `ij` (minus the signed area for edges
that run with the surface on the left), names `h_j := (Uᵀa)_j = ½ (a_k − a_i)` for consecutive
boundary vertices `i, j, k` (Eq. 17), and solves `A b = −Uᵀa`. The right-hand side of that system
is `½ (a_i − a_k)`, predecessor minus successor: this is what `hilbert_transform` returns and what
`hilbertRhs` transcribes. -/

section Hilbert

variable {β : Type*} [Fintype β]

/-- Signed area of the polygon with vertices `(a_j, b_j)`, `j` running along the loop:
`A(a, b) = ½ Σ_j (a_j·b_{next j} − a_{next j}·b_j)`; positive for a counter-clockwise polygon. -/
def loopSignedArea (next : Equiv.Perm β) (a b : β → ℝ) : ℝ :=
  1 / 2 * ∑ j, (a j * b (next j) - a (next j) * b j)

/-- What `hilbert_transform` writes at boundary vertex `j`: `h_j = ½ (a_{prev j} − a_{next j})`.
Rust: `0.5 * (a[previous] - a[next])`. -/
def hilbertRhs {β : Type*} (next : Equiv.Perm β) (a : β → ℝ) (j : β) : ℝ :=
  1 / 2 * (a (next.symm j) - a (next j))

/-- `hilbertRhs` in the Rust indexing. -/
theorem hilbertRhs_finRotate {n : ℕ} (a : ℕ → ℝ) (j : Fin n) :
    hilbertRhs (finRotate n) (fun i => a i) j =
      1 / 2 * (a ((j + n - 1) % n) - a ((j + 1) % n)) := by
  simp only [hilbertRhs, finRotate_symm_val, finRotate_val]

/-- **The signed area is `A(a, b) = Σ_j h_j·b_j`**: the vector `h` of `hilbert_transform` is the
gradient of the area with respect to `b`.

[formula: bff-hilbert-transform] -/
theorem loopSignedArea_eq_sum_hilbertRhs_mul (next : Equiv.Perm β) (a b : β → ℝ) :
    loopSignedArea next a b = ∑ j, hilbertRhs next a j * b j := by
  have h : ∑ j, a j * b (next j) = ∑ j, a (next.symm j) * b j := by
    rw [← Equiv.sum_comp next (fun j => a (next.symm j) * b j)]
    simp
  calc loopSignedArea next a b
      = 1 / 2 * (∑ j, a j * b (next j) - ∑ j, a (next j) * b j) := by
        rw [loopSignedArea, Finset.sum_sub_distrib]
    _ = 1 / 2 * (∑ j, a (next.symm j) * b j - ∑ j, a (next j) * b j) := by rw [h]
    _ = ∑ j, hilbertRhs next a j * b j := by
        rw [← Finset.sum_sub_distrib, Finset.mul_sum]
        exact Finset.sum_congr rfl fun j _ => by simp only [hilbertRhs]; ring

/-- **The right-hand side sums to zero**: the compatibility condition of the Neumann problem.

[formula: bff-hilbert-transform] -/
theorem sum_hilbertRhs (next : Equiv.Perm β) (a : β → ℝ) : ∑ j, hilbertRhs next a j = 0 := by
  simp only [hilbertRhs]
  rw [← Finset.mul_sum, Finset.sum_sub_distrib, Equiv.sum_comp next.symm a, Equiv.sum_comp next a,
    sub_self, mul_zero]

/-- The signed area is antisymmetric: `A(a, b) = −A(b, a)`. -/
theorem loopSignedArea_antisymm (next : Equiv.Perm β) (a b : β → ℝ) :
    loopSignedArea next a b = -loopSignedArea next b a := by
  simp only [loopSignedArea]
  rw [← mul_neg, ← Finset.sum_neg_distrib]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by ring

/-- `A(a, a) = 0`. -/
theorem loopSignedArea_self (next : Equiv.Perm β) (a : β → ℝ) : loopSignedArea next a a = 0 := by
  have h := loopSignedArea_antisymm next a a
  linarith

/-- The unit square `(0,0), (1,0), (1,1), (0,1)`, counter-clockwise, has signed area `+1`. -/
example : loopSignedArea (finRotate 4) ![0, 1, 1, 0] ![0, 0, 1, 1] = 1 := by
  simp [loopSignedArea, Fin.sum_univ_four]
  norm_num

/-! ### The Neumann problem

`V` is the set of vertices of the mesh, `emb : β → V` sends a position on the boundary loop to
its vertex (the array `boundary` of the Rust code), `L` is the cotan Laplacian. -/

variable {V : Type*} [Fintype V]

/-- The part of the conformal energy that depends on `b`: `½ bᵀLb − A(a, b)`, the area being
that of the image of the boundary loop. -/
def conjugateEnergy (next : Equiv.Perm β) (emb : β → V) (L : Matrix V V ℝ) (a b : V → ℝ) : ℝ :=
  1 / 2 * (b ⬝ᵥ (L *ᵥ b)) - loopSignedArea next (a ∘ emb) (b ∘ emb)

/-- The conformal energy of the map `(a, b)`: `E_C(a, b) = ½ (aᵀLa + bᵀLb) − A(a, b)`, the
Dirichlet energy minus the signed area of the image. -/
def loopConformalEnergy (next : Equiv.Perm β) (emb : β → V) (L : Matrix V V ℝ) (a b : V → ℝ) :
    ℝ :=
  1 / 2 * (a ⬝ᵥ (L *ᵥ a) + b ⬝ᵥ (L *ᵥ b)) - loopSignedArea next (a ∘ emb) (b ∘ emb)

/-- The right-hand side as `hilbert_transform` builds it: zero everywhere, then `h_j` written at
the vertex of boundary position `j`. (For a loop that visits no vertex twice. The definition adds
the values written at the same vertex, the Rust code would keep the last one.) -/
def hilbertVector [DecidableEq V] (next : Equiv.Perm β) (emb : β → V) (a : V → ℝ)
    (v : V) : ℝ :=
  ∑ j, if emb j = v then hilbertRhs next (a ∘ emb) j else 0

/-- `E_C(a, b) = ½ aᵀLa + (½ bᵀLb − A(a, b))`. -/
theorem loopConformalEnergy_eq (next : Equiv.Perm β) (emb : β → V) (L : Matrix V V ℝ)
    (a b : V → ℝ) :
    loopConformalEnergy next emb L a b =
      1 / 2 * (a ⬝ᵥ (L *ᵥ a)) + conjugateEnergy next emb L a b := by
  simp only [loopConformalEnergy, conjugateEnergy]
  ring

omit [Fintype V] in
/-- At the vertex of boundary position `j` the vector of `hilbert_transform` is `h_j`, if no
vertex occurs twice in the loop. -/
theorem hilbertVector_emb [DecidableEq V] (next : Equiv.Perm β) {emb : β → V}
    (hemb : Function.Injective emb) (a : V → ℝ) (j : β) :
    hilbertVector next emb a (emb j) = hilbertRhs next (a ∘ emb) j := by
  rw [hilbertVector, Finset.sum_eq_single j]
  · simp
  · intro i _ hij
    rw [ite_eq_right_iff]
    intro h
    exact absurd (hemb h) hij
  · intro h
    exact absurd (Finset.mem_univ j) h

omit [Fintype V] in
/-- At vertices that are not on the loop the vector of `hilbert_transform` is zero. -/
theorem hilbertVector_of_notMem [DecidableEq V] (next : Equiv.Perm β) (emb : β → V) (a : V → ℝ)
    {v : V} (hv : v ∉ Set.range emb) : hilbertVector next emb a v = 0 := by
  refine Finset.sum_eq_zero fun j _ => ?_
  rw [ite_eq_right_iff]
  intro h
  exact absurd ⟨j, h⟩ hv

/-- Pairing the right-hand side with `b` gives the signed area of the image of the boundary:
`rhs · b = A(a, b)`. -/
theorem dotProduct_hilbertRhs {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) (b : V → ℝ) :
    rhs ⬝ᵥ b = loopSignedArea next (a ∘ emb) (b ∘ emb) := by
  rw [loopSignedArea_eq_sum_hilbertRhs_mul, dotProduct]
  symm
  exact Fintype.sum_of_injective emb hemb _ (fun v => rhs v * b v)
    (fun v hv => by simp [hzero v hv]) (fun j => by simp [hrhs j])

/-- The right-hand side sums to zero over all vertices: it is orthogonal to the constants, the
kernel of the Laplacian of a connected mesh. -/
theorem sum_hilbertRhs_vertices {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) : ∑ v, rhs v = 0 := by
  rw [← Fintype.sum_of_injective emb hemb (hilbertRhs next (a ∘ emb)) rhs hzero
    (fun j => (hrhs j).symm)]
  exact sum_hilbertRhs next _

private theorem dotProduct_mulVec_comm_of_isSymm {L : Matrix V V ℝ} (hL : L.IsSymm)
    (x y : V → ℝ) : x ⬝ᵥ (L *ᵥ y) = y ⬝ᵥ (L *ᵥ x) := by
  rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hL.eq, dotProduct_comm]

/-- The energy along a segment: for symmetric `L`,
`E(b') − E(b) = (b' − b)·(L b − rhs) + ½ (b' − b)ᵀ L (b' − b)`. -/
theorem conjugateEnergy_sub {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.IsSymm) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) (b b' : V → ℝ) :
    conjugateEnergy next emb L a b' - conjugateEnergy next emb L a b =
      (b' - b) ⬝ᵥ (L *ᵥ b - rhs) + 1 / 2 * ((b' - b) ⬝ᵥ (L *ᵥ (b' - b))) := by
  have hA := dotProduct_hilbertRhs hemb hrhs hzero
  have hs : b ⬝ᵥ (L *ᵥ b') = b' ⬝ᵥ (L *ᵥ b) := dotProduct_mulVec_comm_of_isSymm hL b b'
  simp only [conjugateEnergy, ← hA, Matrix.mulVec_sub, dotProduct_sub, sub_dotProduct, hs,
    dotProduct_comm rhs]
  ring

/-- **A solution of `L b = rhs` minimizes the conformal energy over `b`** (paper Eq. 17). Here
`L` is symmetric positive semidefinite, `emb` is injective, and `rhs` is `h_j = ½ (a_{prev j} −
a_{next j})` at the vertex of boundary position `j` and zero at the vertices that are not on the
loop. Then `½ bᵀLb − A(a, b) ≤ ½ b'ᵀLb' − A(a, b')` for every `b'`.

[formula: bff-hilbert-transform] -/
theorem hilbertTransform_minimizes {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.IsSymm)
    (hpsd : ∀ x : V → ℝ, 0 ≤ x ⬝ᵥ (L *ᵥ x)) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) {b : V → ℝ} (hb : L *ᵥ b = rhs)
    (b' : V → ℝ) :
    conjugateEnergy next emb L a b ≤ conjugateEnergy next emb L a b' := by
  have h := conjugateEnergy_sub hemb hL hrhs hzero b b'
  rw [hb, sub_self, dotProduct_zero] at h
  have := hpsd (b' - b)
  linarith

/-- **Every minimizer of the conformal energy over `b` solves `L b = rhs`.** Only the symmetry of
`L` is used.

[formula: bff-hilbert-transform] -/
theorem hilbertTransform_of_minimizes {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.IsSymm) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) {b : V → ℝ}
    (hmin : ∀ b', conjugateEnergy next emb L a b ≤ conjugateEnergy next emb L a b') :
    L *ᵥ b = rhs := by
  set g := L *ᵥ b - rhs with hg
  have key : ∀ t : ℝ, 0 ≤ t * (g ⬝ᵥ g) + 1 / 2 * (t ^ 2 * (g ⬝ᵥ (L *ᵥ g))) := by
    intro t
    have h := conjugateEnergy_sub hemb hL hrhs hzero b (b + t • g)
    have := hmin (b + t • g)
    rw [add_sub_cancel_left] at h
    simp only [smul_dotProduct, dotProduct_smul, Matrix.mulVec_smul, smul_eq_mul] at h
    nlinarith
  have hq : g ⬝ᵥ g = 0 := by
    by_contra hne
    have hnn : 0 ≤ g ⬝ᵥ g := Finset.sum_nonneg fun i _ => mul_self_nonneg (g i)
    have hpos : 0 < g ⬝ᵥ g := lt_of_le_of_ne hnn (Ne.symm hne)
    rcases le_or_gt (g ⬝ᵥ (L *ᵥ g)) 0 with hc | hc
    · have := key (-1)
      nlinarith
    · have h := key (-(g ⬝ᵥ g / (g ⬝ᵥ (L *ᵥ g))))
      have e : -(g ⬝ᵥ g / (g ⬝ᵥ (L *ᵥ g))) * (g ⬝ᵥ g) +
          1 / 2 * ((-(g ⬝ᵥ g / (g ⬝ᵥ (L *ᵥ g)))) ^ 2 * (g ⬝ᵥ (L *ᵥ g))) =
          -((g ⬝ᵥ g) ^ 2 / (2 * (g ⬝ᵥ (L *ᵥ g)))) := by
        field_simp
        ring
      rw [e] at h
      have : 0 < (g ⬝ᵥ g) ^ 2 / (2 * (g ⬝ᵥ (L *ᵥ g))) := by positivity
      linarith
  exact sub_eq_zero.1 (dotProduct_self_eq_zero.1 hq)

/-- The minimizers of `b ↦ ½ bᵀLb − A(a, b)` are exactly the solutions of `L b = rhs`. -/
theorem hilbertTransform_minimizes_iff {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.IsSymm)
    (hpsd : ∀ x : V → ℝ, 0 ≤ x ⬝ᵥ (L *ᵥ x)) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) (b : V → ℝ) :
    (∀ b', conjugateEnergy next emb L a b ≤ conjugateEnergy next emb L a b') ↔ L *ᵥ b = rhs :=
  ⟨hilbertTransform_of_minimizes hemb hL hrhs hzero,
    fun hb => hilbertTransform_minimizes hemb hL hpsd hrhs hzero hb⟩

/-- **The minimizers over `b` of the conformal energy `E_C(a, b) = ½ (aᵀLa + bᵀLb) − A(a, b)`
are exactly the solutions of `L b = rhs`**, for `L` symmetric positive semidefinite.

[formula: bff-hilbert-transform] -/
theorem hilbertTransform_minimizes_conformalEnergy_iff {next : Equiv.Perm β} {emb : β → V}
    (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.IsSymm)
    (hpsd : ∀ x : V → ℝ, 0 ≤ x ⬝ᵥ (L *ᵥ x)) {a rhs : V → ℝ}
    (hrhs : ∀ j, rhs (emb j) = hilbertRhs next (a ∘ emb) j)
    (hzero : ∀ v, v ∉ Set.range emb → rhs v = 0) (b : V → ℝ) :
    (∀ b', loopConformalEnergy next emb L a b ≤ loopConformalEnergy next emb L a b') ↔
      L *ᵥ b = rhs := by
  rw [← hilbertTransform_minimizes_iff hemb hL hpsd hrhs hzero b]
  simp only [loopConformalEnergy_eq, add_le_add_iff_left]

/-- The same for the vector that `hilbert_transform` returns and a matrix that is positive
semidefinite in the sense of Mathlib (as the cotan Laplacian is): `b` minimizes the conformal
energy `E_C(a, ·)` if and only if `L b = hilbert_transform(a)`. -/
theorem hilbertVector_minimizes_conformalEnergy_iff [DecidableEq V] (next : Equiv.Perm β)
    {emb : β → V} (hemb : Function.Injective emb) {L : Matrix V V ℝ} (hL : L.PosSemidef)
    (a b : V → ℝ) :
    (∀ b', loopConformalEnergy next emb L a b ≤ loopConformalEnergy next emb L a b') ↔
      L *ᵥ b = hilbertVector next emb a := by
  have hsymm : L.IsSymm := by
    have h := hL.isHermitian
    rwa [Matrix.IsHermitian, Matrix.conjTranspose_eq_transpose_of_trivial] at h
  have hpsd : ∀ x : V → ℝ, 0 ≤ x ⬝ᵥ (L *ᵥ x) := fun x => by
    simpa using hL.dotProduct_mulVec_nonneg x
  exact hilbertTransform_minimizes_conformalEnergy_iff hemb hsymm hpsd
    (hilbertVector_emb next hemb a) (fun v hv => hilbertVector_of_notMem next emb a hv) b

end Hilbert

end DDG
