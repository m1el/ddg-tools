import DdgProofs.Cuts.Recovery

/-!
# Cut flow: the default stages on the cut

`crates/ddg-variational-cuts/SPEC.md` §3.8 – §3.10, D14 – D16: the length force and the rate are
computed on the cut, the rate is carried to the vertices by weighted means, and it is scaled by
the mean of `|∇φ|` over the piece of the cut.

The definitions of this file (`wmean`, `localExtension`, `gradNorm`, `meanGradient`,
`rippleFactor`) are **this package's own**; the two linear systems are `Cuts.fluxMatrix`, built
from the oracle's `fluxDiagEntry`, `fluxOffEntry`, `fluxRhs` (`DdgProofs/Cuts/Recovery.lean`).

* **`cut-length-step-on-cut`** — the two systems of SPEC §3.10,
  `ℓ_p λ⁰_p + Σ κ_σ (λ⁰_p − λ⁰_q) = Θ_p` and
  `ℓ_p r_p + τ Σ (w̄_σ/|pq|) (r_p − r_q) = ℓ_p (s_p + λ⁰_p)`:
  both matrices are symmetric positive definite with nonpositive off-diagonal entries, row sums
  `ℓ_p` and a nonnegative inverse (`fluxMatrix_apply_of_ne_nonpos`, `fluxMatrix_row_sum`,
  `fluxMatrix_posDef`, `fluxMatrix_monotone`), so `λ⁰` and `r` are unique
  (`onCut_force_existsUnique`, `onCut_rate_existsUnique`); `Σ ℓ λ⁰ = Σ Θ` and
  `Σ ℓ r = Σ ℓ (s + λ⁰)` (`onCut_force_sum`, `onCut_rate_sum`, also per piece of the cut); `r`
  lies between the smallest and the largest `s + λ⁰` (`onCut_rate_bounds`); if `s + λ⁰ = 0` then
  `r = 0` for every `τ` (`onCut_rest`); without coupling `r = s + λ⁰` (`onCut_rate_nodal`). The
  per-step factor of a ripple, `1 + τ(ωk − wk²)/(1 + τwk²) = (1 + τωk)/(1 + τwk²)`, lies in
  `(0, 1)` for `k > ω/w`, `τ > 0`, `w > 0` **and `ω ≥ 0`** (`rippleFactor_mem_Ioo`; for `ω < 0`
  it can be negative, `rippleFactor_neg_example`). The turning angles `π − β_p` around a region
  sum to `2πχ − Σ_interior Ω` by discrete Gauss–Bonnet, `2π` for a flat disk
  (`TriMesh.sum_turningAngle`, `TriMesh.onCut_force_sum_flat_disk`): a closed cut on a flat
  surface has `Σ ℓ λ⁰ = 2πw`.
* **`cut-extension-local`** — `σ_i = Σ_p ℓ_p β_ip r_p / Σ_p ℓ_p β_ip` is a convex combination of
  the rates at the cut points of `i` (`localExtension_mem`); the outward rule is one too
  (`outwardExtension_mem`); every vertex that the rings reach has `σ` between the smallest and
  the largest rate (`extension_mem_of_rings`); `Σ_p ℓ_p Σ_i β_ip σ_i = Σ_p ℓ_p r_p` when the two
  barycentric weights of every cut point sum to one (`localExtension_conservation`).
* **`cut-gradient-norm`** — for an affine function on a triangle `a p q` that vanishes at `p`
  and `q` and has the value `φ_a` at `a`, `‖∇φ‖ = |φ_a| · |pq| / (2 · area(a, p, q))`
  (`norm_grad_eq_gradNorm`); the mean `ḡ = Σ|pq| / Σ(|pq|/|∇φ|)` lies between the smallest and
  the largest `|∇φ|` of the piece, and the mean of `ḡ/|∇φ|` over the piece is `1`
  (`meanGradient_mem`, `meanGradient_normalizes`).
-/

noncomputable section

namespace DdgProofs

namespace Cuts

open Matrix
open scoped DdgProofs.Oracle
open Real

/-! ## Weighted means -/

/-- The weighted mean `Σ_{j ∈ s} c_j x_j / Σ_{j ∈ s} c_j`. -/
def wmean {ι : Type*} (s : Finset ι) (c x : ι → ℝ) : ℝ := (∑ j ∈ s, c j * x j) / ∑ j ∈ s, c j

/-- **A weighted mean with nonnegative weights is a convex combination**: it lies between any
lower and upper bound of the values that carry weight. -/
theorem wmean_mem {ι : Type*} {s : Finset ι} {c x : ι → ℝ} (hc : ∀ j ∈ s, 0 ≤ c j)
    (hpos : 0 < ∑ j ∈ s, c j) {lo hi : ℝ} (hlo : ∀ j ∈ s, c j ≠ 0 → lo ≤ x j)
    (hhi : ∀ j ∈ s, c j ≠ 0 → x j ≤ hi) : lo ≤ wmean s c x ∧ wmean s c x ≤ hi := by
  unfold wmean
  rw [le_div_iff₀ hpos, div_le_iff₀ hpos, Finset.mul_sum, Finset.mul_sum]
  constructor
  · refine Finset.sum_le_sum fun j hj => ?_
    by_cases h0 : c j = 0
    · rw [h0]
      simp
    · rw [mul_comm]
      exact mul_le_mul_of_nonneg_left (hlo j hj h0) (hc j hj)
  · refine Finset.sum_le_sum fun j hj => ?_
    by_cases h0 : c j = 0
    · rw [h0]
      simp
    · rw [mul_comm hi]
      exact mul_le_mul_of_nonneg_left (hhi j hj h0) (hc j hj)

/-- The weighted mean of a constant is that constant. -/
theorem wmean_const {ι : Type*} {s : Finset ι} {c : ι → ℝ} (hne : ∑ j ∈ s, c j ≠ 0) (a : ℝ) :
    wmean s c (fun _ => a) = a := by
  unfold wmean
  rw [← Finset.sum_mul]
  field_simp

/-- The weighted mean is linear in the values. -/
theorem wmean_smul {ι : Type*} (s : Finset ι) (c x : ι → ℝ) (a : ℝ) :
    wmean s c (fun j => a * x j) = a * wmean s c x := by
  unfold wmean
  rw [← mul_div_assoc, Finset.mul_sum]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by ring

theorem wmean_add {ι : Type*} (s : Finset ι) (c x y : ι → ℝ) :
    wmean s c (fun j => x j + y j) = wmean s c x + wmean s c y := by
  unfold wmean
  rw [← add_div, ← Finset.sum_add_distrib]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by ring

/-! ## Extension of the rate from the cut (SPEC §3.8, D15) -/

section Extension

variable {V P : Type*} [Fintype V] [Fintype P]

/-- The rate at a cut-adjacent vertex `i` from the rates `r_p` at the cut points:
`σ_i = Σ_p ℓ_p β_ip r_p / Σ_p ℓ_p β_ip`, with `ℓ_p` the dual length of `p` on the cut and `β_ip`
the barycentric weight of `i` at `p` (`1 − t` or `t` on the edge of `p`, `0` if `p` is not on an
edge at `i`). -/
def localExtension (ℓ : P → ℝ) (β : V → P → ℝ) (r : P → ℝ) (i : V) : ℝ :=
  wmean Finset.univ (fun p => ℓ p * β i p) r

omit [Fintype V] in
/-- **The value at a cut-adjacent vertex is a convex combination of the rates at its cut
points**: it lies between the smallest and the largest of them.

[formula: cut-extension-local] -/
theorem localExtension_mem {ℓ : P → ℝ} {β : V → P → ℝ} (hℓ : ∀ p, 0 ≤ ℓ p)
    (hβ : ∀ i p, 0 ≤ β i p) {r : P → ℝ} {i : V} (hi : 0 < ∑ p, ℓ p * β i p) {lo hi' : ℝ}
    (hlo : ∀ p, ℓ p * β i p ≠ 0 → lo ≤ r p) (hhi : ∀ p, ℓ p * β i p ≠ 0 → r p ≤ hi') :
    lo ≤ localExtension ℓ β r i ∧ localExtension ℓ β r i ≤ hi' :=
  wmean_mem (fun p _ => mul_nonneg (hℓ p) (hβ i p)) hi (fun p _ => hlo p) (fun p _ => hhi p)

/-- **The mean rate of the cut is kept.** The explicit step interpolates the vertex rates back
to the cut points with the same weights, `Σ_i β_ip σ_i`; if the weights of every cut point sum to
one (`(1 − t) + t`), then `Σ_p ℓ_p Σ_i β_ip σ_i = Σ_p ℓ_p r_p`.

[formula: cut-extension-local] -/
theorem localExtension_conservation {ℓ : P → ℝ} {β : V → P → ℝ} (hℓ : ∀ p, 0 ≤ ℓ p)
    (hβ : ∀ i p, 0 ≤ β i p) (hsum : ∀ p, ∑ i, β i p = 1) (r : P → ℝ) :
    ∑ p, ℓ p * ∑ i, β i p * localExtension ℓ β r i = ∑ p, ℓ p * r p := by
  have hrow : ∀ i, localExtension ℓ β r i * ∑ p, ℓ p * β i p = ∑ p, ℓ p * β i p * r p := by
    intro i
    unfold localExtension wmean
    by_cases h0 : ∑ p, ℓ p * β i p = 0
    · have hz : ∀ p, ℓ p * β i p = 0 := fun p =>
        (Finset.sum_eq_zero_iff_of_nonneg fun p _ => mul_nonneg (hℓ p) (hβ i p)).mp h0 p
          (Finset.mem_univ p)
      rw [h0, mul_zero]
      exact (Finset.sum_eq_zero fun p _ => by rw [hz p, zero_mul]).symm
    · rw [div_mul_cancel₀ _ h0]
  calc ∑ p, ℓ p * ∑ i, β i p * localExtension ℓ β r i
      = ∑ i, localExtension ℓ β r i * ∑ p, ℓ p * β i p := by
        simp only [Finset.mul_sum]
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun p _ => by ring
    _ = ∑ i, ∑ p, ℓ p * β i p * r p := Finset.sum_congr rfl fun i _ => hrow i
    _ = ∑ p, ℓ p * r p := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun p _ => ?_
        have : ∑ i, ℓ p * β i p * r p = ℓ p * r p * ∑ i, β i p := by
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun i _ => by ring
        rw [this, hsum p, mul_one]

omit [Fintype V] [Fintype P] in
/-- **The outward rule is a convex combination too**:
`σ_i = Σ_j (σ_j/ℓ_ij) / Σ_j (1/ℓ_ij)` over a nonempty set `N` of neighbours with positive edge
lengths lies between the smallest and the largest `σ_j`, `j ∈ N`.

[formula: cut-extension-local] -/
theorem outwardExtension_mem {N : Finset V} (hN : N.Nonempty) {len σ : V → ℝ}
    (hlen : ∀ j ∈ N, 0 < len j) {lo hi : ℝ} (hlo : ∀ j ∈ N, lo ≤ σ j)
    (hhi : ∀ j ∈ N, σ j ≤ hi) :
    lo ≤ (∑ j ∈ N, σ j / len j) / ∑ j ∈ N, 1 / len j ∧
      (∑ j ∈ N, σ j / len j) / ∑ j ∈ N, 1 / len j ≤ hi := by
  have e : (∑ j ∈ N, σ j / len j) / ∑ j ∈ N, 1 / len j = wmean N (fun j => 1 / len j) σ := by
    unfold wmean
    congr 1
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [e]
  exact wmean_mem (fun j hj => (one_div_pos.mpr (hlen j hj)).le)
    (Finset.sum_pos (fun j hj => one_div_pos.mpr (hlen j hj)) hN)
    (fun j hj _ => hlo j hj) (fun j hj _ => hhi j hj)

omit [Fintype P] in
/-- **Every vertex that the rings reach gets a value between the smallest and the largest
rate.** `ring i` is the ring of the vertex; ring `0` (cut-adjacent) has values in `[lo, hi]`
(`localExtension_mem`), and a vertex of ring `n + 1` takes a weighted mean, with nonnegative
weights, of vertices of ring `n`. (A vertex that no ring reaches gets `0`, which need not lie in
that range.)

[formula: cut-extension-local] -/
theorem extension_mem_of_rings (σ : V → ℝ) (ring : V → ℕ) {lo hi : ℝ}
    (h0 : ∀ i, ring i = 0 → lo ≤ σ i ∧ σ i ≤ hi)
    (hstep : ∀ i n, ring i = n + 1 → ∃ c : V → ℝ, (∀ j, 0 ≤ c j) ∧ 0 < ∑ j, c j ∧
      (∀ j, c j ≠ 0 → ring j = n) ∧ σ i = wmean Finset.univ c σ) (i : V) :
    lo ≤ σ i ∧ σ i ≤ hi := by
  have aux : ∀ n i, ring i = n → lo ≤ σ i ∧ σ i ≤ hi := by
    intro n
    induction n with
    | zero => exact h0
    | succ n ih =>
      intro i hi
      obtain ⟨c, hc0, hcpos, hcring, hσ⟩ := hstep i n hi
      rw [hσ]
      exact wmean_mem (fun j _ => hc0 j) hcpos (fun j _ hj => (ih j (hcring j hj)).1)
        (fun j _ hj => (ih j (hcring j hj)).2)
  exact aux (ring i) i rfl

end Extension

/-! ## The gradient of the level-set function in a crossed face (SPEC §3.9, D16) -/

/-- `|∇φ|` of the linear interpolant in a face whose lone vertex `a` is cut off by the segment
`pq`: `|φ_a| · |pq| / (2 · area(a, p, q))`. -/
def gradNorm (φa len area : ℝ) : ℝ := |φa| * len / (2 * area)

/-- **The gradient of a linear function that vanishes on the cut segment.** Triangle `a p q`
with `e₁ = p − a`, `e₂ = q − a` (nondegenerate); `g = α e₁ + β e₂` the gradient of the affine
function `φ(x) = φ_a + ⟨g, x − a⟩` in the plane of the triangle; `φ(p) = φ(q) = 0`, i.e.
`⟨g, e₁⟩ = ⟨g, e₂⟩ = −φ_a`. Then `‖g‖ = |φ_a| · |pq| / (2 · area(a, p, q))`, the area being
Heron's of the three side lengths.

[formula: cut-gradient-norm] -/
theorem norm_grad_eq_gradNorm {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (e₁ e₂ : E) (α β φa : ℝ) (h : IsTriangle ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖)
    (h1 : inner ℝ (α • e₁ + β • e₂) e₁ = -φa) (h2 : inner ℝ (α • e₁ + β • e₂) e₂ = -φa) :
    ‖α • e₁ + β • e₂‖ = gradNorm φa ‖e₂ - e₁‖ (triArea ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖) := by
  have hT := sixteen_mul_triArea_sq h
  have hTpos := triArea_pos h
  have hc2 : inner ℝ e₁ e₁ = ‖e₁‖ ^ 2 := real_inner_self_eq_norm_sq e₁
  have hb2 : inner ℝ e₂ e₂ = ‖e₂‖ ^ 2 := real_inner_self_eq_norm_sq e₂
  have hsymm : inner ℝ e₂ e₁ = inner ℝ e₁ e₂ := real_inner_comm _ _
  have ha2 : ‖e₂ - e₁‖ ^ 2 = ‖e₂‖ ^ 2 - 2 * inner ℝ e₁ e₂ + ‖e₁‖ ^ 2 := by
    rw [norm_sub_sq_real, hsymm]
  have hg : ‖α • e₁ + β • e₂‖ ^ 2 =
      α ^ 2 * ‖e₁‖ ^ 2 + 2 * α * β * inner ℝ e₁ e₂ + β ^ 2 * ‖e₂‖ ^ 2 := by
    rw [← real_inner_self_eq_norm_sq]
    simp only [inner_add_left, inner_add_right, inner_smul_left, inner_smul_right, hc2, hb2,
      hsymm, RCLike.conj_to_real]
    ring
  have hf1 : α * ‖e₁‖ ^ 2 + β * inner ℝ e₁ e₂ = -φa := by
    rw [← h1]
    simp only [inner_add_left, inner_smul_left, hc2, hsymm, RCLike.conj_to_real]
  have hf2 : α * inner ℝ e₁ e₂ + β * ‖e₂‖ ^ 2 = -φa := by
    rw [← h2]
    simp only [inner_add_left, inner_smul_left, hb2, RCLike.conj_to_real]
  have hn := norm_nonneg (α • e₁ + β • e₂)
  have hl := norm_nonneg (e₂ - e₁)
  generalize ‖α • e₁ + β • e₂‖ = n at *
  generalize inner ℝ e₁ e₂ = F at *
  generalize ‖e₂ - e₁‖ = a at *
  generalize ‖e₂‖ = b at *
  generalize ‖e₁‖ = c at *
  generalize triArea a b c = T at *
  have h4 : a ^ 4 = (b ^ 2 - 2 * F + c ^ 2) ^ 2 := by rw [← ha2]; ring
  have hT' : 4 * T ^ 2 = b ^ 2 * c ^ 2 - F ^ 2 := by
    linear_combination (1 / 4) * hT + (1 / 2) * (b ^ 2 + c ^ 2) * ha2 - (1 / 4) * h4
  have key : n ^ 2 * (b ^ 2 * c ^ 2 - F ^ 2) = φa ^ 2 * (b ^ 2 - 2 * F + c ^ 2) := by
    rw [hg]
    linear_combination (b ^ 2 * (α * c ^ 2 + β * F - φa) - 2 * F * (α * F + β * b ^ 2)) * hf1 +
      (2 * F * φa + c ^ 2 * (α * F + β * b ^ 2 - φa)) * hf2
  have hsq : (n * (2 * T)) ^ 2 = (|φa| * a) ^ 2 := by
    have e1 : (|φa| * a) ^ 2 = φa ^ 2 * a ^ 2 := by rw [mul_pow, sq_abs]
    rw [e1]
    linear_combination key + n ^ 2 * hT' - φa ^ 2 * ha2
  have heq : n * (2 * T) = |φa| * a := by
    have h1' : 0 ≤ n * (2 * T) := by positivity
    have h2' : 0 ≤ |φa| * a := by positivity
    calc n * (2 * T) = √((n * (2 * T)) ^ 2) := (Real.sqrt_sq h1').symm
      _ = √((|φa| * a) ^ 2) := by rw [hsq]
      _ = |φa| * a := Real.sqrt_sq h2'
  unfold gradNorm
  rw [eq_div_iff (by positivity)]
  exact heq

/-- `gradNorm` is the quotient of `|φ_a|` and the height of `a` over `pq` (`2·area/|pq|`). -/
theorem gradNorm_eq_div_height (φa : ℝ) {len area : ℝ} (hlen : len ≠ 0) (harea : area ≠ 0) :
    gradNorm φa len area = |φa| / (2 * area / len) := by
  unfold gradNorm
  field_simp

/-- `gradNorm` is positive when `φ_a ≠ 0` and the segment and the area are positive. -/
theorem gradNorm_pos {φa len area : ℝ} (hφ : φa ≠ 0) (hlen : 0 < len) (harea : 0 < area) :
    0 < gradNorm φa len area := by
  unfold gradNorm
  have := abs_pos.mpr hφ
  positivity

/-- The mean of `|∇φ|` over a piece of the cut (SPEC §3.9): the length-weighted harmonic mean
`ḡ = Σ |pq| / Σ (|pq| / |∇φ|)` over its segments. -/
def meanGradient {ι : Type*} (s : Finset ι) (L g : ι → ℝ) : ℝ :=
  (∑ σ ∈ s, L σ) / ∑ σ ∈ s, L σ / g σ

/-- **The mean `ḡ` lies between the smallest and the largest `|∇φ|` of the piece.**

[formula: cut-gradient-norm] -/
theorem meanGradient_mem {ι : Type*} {s : Finset ι} {L g : ι → ℝ} (hL : ∀ σ ∈ s, 0 ≤ L σ)
    (hpos : 0 < ∑ σ ∈ s, L σ) {lo hi : ℝ} (hlo0 : 0 < lo) (hlo : ∀ σ ∈ s, lo ≤ g σ)
    (hhi : ∀ σ ∈ s, g σ ≤ hi) : lo ≤ meanGradient s L g ∧ meanGradient s L g ≤ hi := by
  have hg0 : ∀ σ ∈ s, 0 < g σ := fun σ hσ => hlo0.trans_le (hlo σ hσ)
  have h1 : lo * ∑ σ ∈ s, L σ / g σ ≤ ∑ σ ∈ s, L σ := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun σ hσ => ?_
    have hgσ := hg0 σ hσ
    rw [← mul_div_assoc, div_le_iff₀ hgσ, mul_comm]
    exact mul_le_mul_of_nonneg_left (hlo σ hσ) (hL σ hσ)
  have h2 : ∑ σ ∈ s, L σ ≤ hi * ∑ σ ∈ s, L σ / g σ := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum fun σ hσ => ?_
    have hgσ := hg0 σ hσ
    rw [← mul_div_assoc, le_div_iff₀ hgσ, mul_comm hi]
    exact mul_le_mul_of_nonneg_left (hhi σ hσ) (hL σ hσ)
  have hhi0 : 0 < hi := by
    obtain ⟨σ, hσ, -⟩ : ∃ σ ∈ s, 0 < L σ := by
      by_contra hcon
      push Not at hcon
      have : ∑ σ ∈ s, L σ ≤ 0 := Finset.sum_nonpos hcon
      linarith
    exact (hg0 σ hσ).trans_le (hhi σ hσ)
  have hT : 0 < ∑ σ ∈ s, L σ / g σ := by
    by_contra hcon
    push Not at hcon
    have := mul_nonpos_of_nonneg_of_nonpos hhi0.le hcon
    linarith
  unfold meanGradient
  rw [le_div_iff₀ hT, div_le_iff₀ hT]
  exact ⟨h1, h2⟩

/-- **Why the harmonic mean**: a rate `r·ḡ` of the function moves a segment with gradient `g_σ`
by `r ḡ/g_σ` per unit time; the length-weighted mean of `ḡ/g_σ` over the piece is `1`, so the
mean normal speed of the piece is `r`.

[formula: cut-gradient-norm] -/
theorem meanGradient_normalizes {ι : Type*} {s : Finset ι} {L g : ι → ℝ}
    (hne : ∑ σ ∈ s, L σ / g σ ≠ 0) :
    ∑ σ ∈ s, L σ * (meanGradient s L g / g σ) = ∑ σ ∈ s, L σ := by
  have e : ∀ σ ∈ s, L σ * (meanGradient s L g / g σ) = meanGradient s L g * (L σ / g σ) :=
    fun σ _ => by ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum]
  unfold meanGradient
  rw [div_mul_cancel₀ _ hne]

/-- If `|∇φ| = 1` on every segment of the piece (a signed distance function), `ḡ = 1`: the
paper's rule `SpeedScaling::Unit`. -/
theorem meanGradient_one {ι : Type*} {s : Finset ι} {L : ι → ℝ} (hpos : ∑ σ ∈ s, L σ ≠ 0) :
    meanGradient s L (fun _ => 1) = 1 := by
  unfold meanGradient
  simp only [div_one]
  exact div_self hpos

/-! ## The two systems on the cut (SPEC §3.10, D14) -/

section OnCut

variable {C : Type*} [Fintype C] [DecidableEq C]

/-- Off-diagonal entries of the flux matrix: `−k(p, q)`. -/
theorem fluxMatrix_apply_of_ne (ℓ : C → ℝ) (k : C → C → ℝ) {p q : C} (h : p ≠ q) :
    fluxMatrix ℓ k p q = -k p q := by
  simp [fluxMatrix, h, fluxOffEntry_eq]

/-- **The off-diagonal entries are nonpositive** (a Z-matrix).

[formula: cut-length-step-on-cut] -/
theorem fluxMatrix_apply_of_ne_nonpos (ℓ : C → ℝ) {k : C → C → ℝ} (hk0 : ∀ p q, 0 ≤ k p q)
    {p q : C} (h : p ≠ q) : fluxMatrix ℓ k p q ≤ 0 := by
  rw [fluxMatrix_apply_of_ne ℓ k h]
  exact neg_nonpos.mpr (hk0 p q)

/-- **The row sums are the dual lengths `ℓ_p`.**

[formula: cut-length-step-on-cut] -/
theorem fluxMatrix_row_sum (ℓ : C → ℝ) (k : C → C → ℝ) (p : C) :
    ∑ q, fluxMatrix ℓ k p q = ℓ p := by
  have h := fluxMatrix_mulVec_apply ℓ k (fun _ => 1) p
  simp only [Matrix.mulVec, dotProduct, mul_one, sub_self, mul_zero, Finset.sum_const_zero,
    add_zero] at h
  exact h

/-- **The inverse is nonnegative** (with the two facts above and positive definiteness: an
M-matrix): if `fluxMatrix · d ≥ 0` entrywise then `d ≥ 0`.

[formula: cut-length-step-on-cut] -/
theorem fluxMatrix_monotone {ℓ : C → ℝ} (hℓ : ∀ p, 0 < ℓ p) {k : C → C → ℝ}
    (hk0 : ∀ p q, 0 ≤ k p q) {d : C → ℝ} (h : ∀ p, 0 ≤ (fluxMatrix ℓ k *ᵥ d) p) (p : C) :
    0 ≤ d p :=
  fluxRecovery_ge_piece hℓ hk0 (v := fluxMatrix ℓ k *ᵥ d) (funext fun _ => rfl) Finset.univ
    (fun _ _ q hq => absurd (Finset.mem_univ q) hq)
    (fun q _ => div_nonneg (h q) (hℓ q).le) p (Finset.mem_univ p)

/-- **The force per unit length `λ⁰` exists and is unique**:
`ℓ_p λ⁰_p + Σ_σ κ_σ (λ⁰_p − λ⁰_q) = Θ_p` (SPEC §3.10, the system of the flux recovery with the
right-hand side `Θ`).

[formula: cut-length-step-on-cut] -/
theorem onCut_force_existsUnique {ℓ : C → ℝ} (hℓ : ∀ p, 0 < ℓ p) {κ : C → C → ℝ}
    (hκ : ∀ p q, κ p q = κ q p) (hκ0 : ∀ p q, 0 ≤ κ p q) (Θ : C → ℝ) :
    ∃! lam : C → ℝ, fluxMatrix ℓ κ *ᵥ lam = fun p => Ddg.fluxRhs (Θ p) :=
  fluxRecovery_existsUnique hℓ hκ hκ0 Θ

/-- **The rate `r` exists and is unique**:
`ℓ_p r_p + τ Σ_σ m_σ (r_p − r_q) = ℓ_p (s_p + λ⁰_p)` with `m_σ = w̄_σ/|pq| ≥ 0`, `τ ≥ 0`.

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_existsUnique {ℓ : C → ℝ} (hℓ : ∀ p, 0 < ℓ p) {m : C → C → ℝ}
    (hm : ∀ p q, m p q = m q p) (hm0 : ∀ p q, 0 ≤ m p q) {τ : ℝ} (hτ : 0 ≤ τ) (s lam : C → ℝ) :
    ∃! r : C → ℝ, fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r =
      fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)) :=
  fluxRecovery_existsUnique hℓ (fun p q => by rw [hm p q]) (fun p q => mul_nonneg hτ (hm0 p q))
    _

/-- **The total length force is kept**: `Σ_p ℓ_p λ⁰_p = Σ_p Θ_p` (on every piece of the cut:
`fluxRecovery_sum_piece`).

[formula: cut-length-step-on-cut] -/
theorem onCut_force_sum {ℓ Θ lam : C → ℝ} {κ : C → C → ℝ} (hκ : ∀ p q, κ p q = κ q p)
    (h : fluxMatrix ℓ κ *ᵥ lam = fun p => Ddg.fluxRhs (Θ p)) :
    ∑ p, ℓ p * lam p = ∑ p, Θ p :=
  fluxRecovery_sum hκ h

/-- **The mean rate is kept**: `Σ_p ℓ_p r_p = Σ_p ℓ_p (s_p + λ⁰_p)`, on every piece `S` of the
cut (a set of cut points that no coupling leaves).

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_sum_piece {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ}
    (hm : ∀ p q, m p q = m q p)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (S : Finset C) (hS : ∀ p ∈ S, ∀ q, q ∉ S → m p q = 0) :
    ∑ p ∈ S, ℓ p * r p = ∑ p ∈ S, ℓ p * (s p + lam p) :=
  fluxRecovery_sum_piece (fun p q => by rw [hm p q]) h S
    (fun p hp q hq => by rw [hS p hp q hq, mul_zero])

/-- **The mean rate is kept** (whole cut).

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_sum {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ} (hm : ∀ p q, m p q = m q p)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r =
      fun p => Ddg.fluxRhs (ℓ p * (s p + lam p))) :
    ∑ p, ℓ p * r p = ∑ p, ℓ p * (s p + lam p) :=
  onCut_rate_sum_piece hm h Finset.univ fun _ _ q hq => absurd (Finset.mem_univ q) hq

/-- **The rate lies between the smallest and the largest `s + λ⁰`**, on every piece of the cut.

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_bounds_piece {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ}
    (hℓ : ∀ p, 0 < ℓ p) (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (S : Finset C) (hS : ∀ p ∈ S, ∀ q, q ∉ S → m p q = 0) {lo hi : ℝ}
    (hlo : ∀ p ∈ S, lo ≤ s p + lam p) (hhi : ∀ p ∈ S, s p + lam p ≤ hi) :
    ∀ p ∈ S, lo ≤ r p ∧ r p ≤ hi := by
  have e : ∀ q, ℓ q * (s q + lam q) / ℓ q = s q + lam q := fun q => by
    have := (hℓ q).ne'
    field_simp
  have hk0 : ∀ p q, 0 ≤ τ * m p q := fun p q => mul_nonneg hτ (hm0 p q)
  have hS' : ∀ p ∈ S, ∀ q, q ∉ S → τ * m p q = 0 := fun p hp q hq => by
    rw [hS p hp q hq, mul_zero]
  intro p hp
  exact ⟨fluxRecovery_ge_piece hℓ hk0 h S hS' (fun q hq => by
      have := hlo q hq
      rw [← e q] at this
      exact this) p hp,
    fluxRecovery_le_piece hℓ hk0 h S hS' (fun q hq => by
      have := hhi q hq
      rw [← e q] at this
      exact this) p hp⟩

/-- **The rate lies between the smallest and the largest `s + λ⁰`** (whole cut).

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_bounds {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    {lo hi : ℝ} (hlo : ∀ p, lo ≤ s p + lam p) (hhi : ∀ p, s p + lam p ≤ hi) (p : C) :
    lo ≤ r p ∧ r p ≤ hi :=
  onCut_rate_bounds_piece hℓ hm0 hτ h Finset.univ
    (fun _ _ q hq => absurd (Finset.mem_univ q) hq) (fun q _ => hlo q) (fun q _ => hhi q) p
    (Finset.mem_univ p)

/-- **A cut in equilibrium does not move, whatever the step**: if `s + λ⁰ = 0` at every cut
point then `r = 0`, for every `τ ≥ 0`. So the rest states of the flow do not depend on `τ`.

[formula: cut-length-step-on-cut] -/
theorem onCut_rest {ℓ s lam r : C → ℝ} {m : C → C → ℝ} {τ : ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hm0 : ∀ p q, 0 ≤ m p q) (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => τ * m p q) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (h0 : ∀ p, s p + lam p = 0) (p : C) : r p = 0 := by
  have := onCut_rate_bounds hℓ hm0 hτ h (lo := 0) (hi := 0) (fun q => (h0 q).ge)
    (fun q => (h0 q).le) p
  exact le_antisymm this.2 this.1

/-- **Without coupling the rate is `s + λ⁰`** (`τ = 0`, or `w ≡ 0` so that every `m_σ = 0`). -/
theorem onCut_rate_nodal {ℓ s lam r : C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (h : fluxMatrix ℓ (fun _ _ => 0) *ᵥ r = fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)))
    (p : C) : r p = s p + lam p := by
  rw [fluxRecovery_nodal hℓ h p]
  unfold Ddg.normalDerivative
  have := (hℓ p).ne'
  field_simp

/-- With a vanishing length force (`Θ = 0`, e.g. `w ≡ 0`) the force per unit length is `0`. -/
theorem onCut_force_zero {ℓ lam : C → ℝ} {κ : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hκ0 : ∀ p q, 0 ≤ κ p q) (h : fluxMatrix ℓ κ *ᵥ lam = fun _ => Ddg.fluxRhs 0) (p : C) :
    lam p = 0 := by
  have := fluxRecovery_bounds hℓ hκ0 h (lo := 0) (hi := 0) (fun q => by simp)
    (fun q => by simp) p
  exact le_antisymm this.2 this.1

/-- **Opposite loads give opposite solutions**: the two sides of a curve see the same chain
(`ℓ`, couplings) and opposite right-hand sides, so their solutions are opposite. -/
theorem fluxRecovery_neg {ℓ v d d' : C → ℝ} {k : C → C → ℝ} (hℓ : ∀ p, 0 < ℓ p)
    (hk : ∀ p q, k p q = k q p) (hk0 : ∀ p q, 0 ≤ k p q)
    (h : fluxMatrix ℓ k *ᵥ d = fun p => Ddg.fluxRhs (v p))
    (h' : fluxMatrix ℓ k *ᵥ d' = fun p => Ddg.fluxRhs (-v p)) : d' = -d := by
  have h2 := fluxRecovery_smul h (-1)
  have e : (fun p => Ddg.fluxRhs (-1 * v p)) = fun p => Ddg.fluxRhs (-v p) := by
    funext p
    rw [neg_one_mul]
  rw [e] at h2
  have := (fluxRecovery_existsUnique hℓ hk hk0 (fun p => -v p)).unique h' h2
  rw [this]
  simp

end OnCut

/-! ## The per-step factor of a ripple (SPEC §3.10, §1) -/

/-- Growth rate of a ripple of wavenumber `k` of a straight cut (SPEC §1):
`ω k − w k²`. -/
def rippleRate (ω w k : ℝ) : ℝ := ω * k - w * k ^ 2

/-- The factor by which one step of the default stages multiplies a ripple of wavenumber `k`
(SPEC §3.10): `1 + τ (ωk − wk²)/(1 + τ w k²)`. -/
def rippleFactor (τ ω w k : ℝ) : ℝ := 1 + τ * rippleRate ω w k / (1 + τ * w * k ^ 2)

/-- Where the factor comes from: for the mode `k` the rate equation of SPEC §3.10 is
`r (1 + τ w k²) = (ω k − w k²) δ`, and the step `δ ↦ δ + τ r` multiplies `δ` by the factor. -/
theorem ripple_step {τ ω w k r δ : ℝ} (h0 : 1 + τ * w * k ^ 2 ≠ 0)
    (h : r * (1 + τ * w * k ^ 2) = rippleRate ω w k * δ) :
    δ + τ * r = rippleFactor τ ω w k * δ := by
  unfold rippleFactor
  have hr : r = rippleRate ω w k * δ / (1 + τ * w * k ^ 2) := by
    rw [eq_div_iff h0, h]
  rw [hr]
  field_simp

/-- The factor is `(1 + τωk)/(1 + τwk²)`. -/
theorem rippleFactor_eq {τ ω w k : ℝ} (h0 : 1 + τ * w * k ^ 2 ≠ 0) :
    rippleFactor τ ω w k = (1 + τ * ω * k) / (1 + τ * w * k ^ 2) := by
  unfold rippleFactor rippleRate
  rw [eq_div_iff h0, add_mul, div_mul_cancel₀ _ h0]
  ring

/-- **Ripples shorter than `2πw/ω` are damped for every step size**: for `τ > 0`, `w > 0`,
`ω ≥ 0` and `k > ω/w` the factor lies in `(0, 1)`.

[formula: cut-length-step-on-cut] -/
theorem rippleFactor_mem_Ioo {τ ω w k : ℝ} (hτ : 0 < τ) (hw : 0 < w) (hω : 0 ≤ ω)
    (hk : ω / w < k) : 0 < rippleFactor τ ω w k ∧ rippleFactor τ ω w k < 1 := by
  have hk0 : 0 < k := lt_of_le_of_lt (div_nonneg hω hw.le) hk
  have hωk : ω < w * k := by
    rw [div_lt_iff₀ hw] at hk
    linarith
  have hden : 0 < 1 + τ * w * k ^ 2 := by positivity
  rw [rippleFactor_eq hden.ne', lt_div_iff₀ hden, div_lt_iff₀ hden]
  constructor
  · have : 0 ≤ τ * ω * k := by positivity
    linarith
  · have : τ * k * ω < τ * k * (w * k) := mul_lt_mul_of_pos_left hωk (by positivity)
    nlinarith

/-- Ripples longer than `2πw/ω` are amplified by every step: for `0 < k < ω/w` the factor
exceeds `1`. -/
theorem one_lt_rippleFactor {τ ω w k : ℝ} (hτ : 0 < τ) (hw : 0 < w) (hk0 : 0 < k)
    (hk : k < ω / w) : 1 < rippleFactor τ ω w k := by
  have hωk : w * k < ω := by
    rw [lt_div_iff₀ hw] at hk
    linarith
  have hden : 0 < 1 + τ * w * k ^ 2 := by positivity
  rw [rippleFactor_eq hden.ne', lt_div_iff₀ hden]
  have : τ * k * (w * k) < τ * k * ω := mul_lt_mul_of_pos_left hωk (by positivity)
  nlinarith

/-- **The hypothesis `ω ≥ 0` is needed**: with `ω = −10`, `w = 1`, `k = 1 > ω/w` and `τ = 1` the
factor is negative. (`ω = 2(Q_0 + Q_1)` is nonnegative for the Dirichlet energy; SPEC states the
range `(0, 1)` without the hypothesis.) -/
theorem rippleFactor_neg_example : rippleFactor 1 (-10) 1 1 < 0 := by
  unfold rippleFactor rippleRate
  norm_num

/-! ## The turning of a closed cut -/

section Turning

open TriMesh

variable {V F : Type*} [DecidableEq V] [Fintype V] [Fintype F] (M : TriMesh V F)

/-- **The turning angles around a region sum to `2πχ` minus its total curvature.** For the cut
mesh of a region (a valid manifold triangle mesh), `π − β_v` with `β_v` the sum of the corner
angles at the boundary vertex `v` is the turning angle of SPEC §3.10, and by discrete
Gauss–Bonnet `Σ_{v ∈ ∂} (π − β_v) = 2πχ − Σ_{v interior} Ω_v`.

[formula: cut-length-step-on-cut] -/
theorem _root_.DdgProofs.TriMesh.sum_turningAngle (hM : M.Valid) (hman : M.Manifold) :
    ∑ v ∈ M.boundaryVertices, (π - M.angleSum v) =
      2 * π * (M.eulerChar : ℝ) - ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v := by
  have h := M.gauss_bonnet hM hman
  simp only [TriMesh.geodesicCurvature] at h
  linarith

/-- **A closed cut on a flat surface**: if the region is a flat disk (no angle defect at its
interior vertices, Euler characteristic `1`), its boundary turns by `2π` in total. -/
theorem _root_.DdgProofs.TriMesh.sum_turningAngle_flat_disk (hM : M.Valid) (hman : M.Manifold)
    (hflat : ∀ v ∈ M.boundaryVerticesᶜ, M.angleDefect v = 0) (hχ : M.eulerChar = 1) :
    ∑ v ∈ M.boundaryVertices, (π - M.angleSum v) = 2 * π := by
  rw [M.sum_turningAngle hM hman, Finset.sum_eq_zero hflat, hχ]
  push_cast
  ring

/-- **The total length force of a closed cut with a constant weight**: with
`Θ_p = w (π − β_p)` at the boundary vertices of the cut mesh of the region,
`Σ_p ℓ_p λ⁰_p = w (2πχ − Σ_interior Ω)`.

[formula: cut-length-step-on-cut] -/
theorem _root_.DdgProofs.TriMesh.onCut_force_sum_const (hM : M.Valid) (hman : M.Manifold)
    (w : ℝ) {ℓ lam : {v // v ∈ M.boundaryVertices} → ℝ}
    {κ : {v // v ∈ M.boundaryVertices} → {v // v ∈ M.boundaryVertices} → ℝ}
    (hκ : ∀ p q, κ p q = κ q p)
    (h : fluxMatrix ℓ κ *ᵥ lam = fun p => Ddg.fluxRhs (w * (π - M.angleSum p.1))) :
    ∑ p, ℓ p * lam p =
      w * (2 * π * (M.eulerChar : ℝ) - ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v) := by
  rw [onCut_force_sum hκ h, ← Finset.mul_sum, ← M.sum_turningAngle hM hman]
  congr 1
  exact Finset.sum_coe_sort M.boundaryVertices fun v => π - M.angleSum v

/-- **A closed cut on a flat surface has the total length force `2πw`** (so it loses area at
the rate `2πw`, SPEC §3.10).

[formula: cut-length-step-on-cut] -/
theorem _root_.DdgProofs.TriMesh.onCut_force_sum_flat_disk (hM : M.Valid) (hman : M.Manifold)
    (hflat : ∀ v ∈ M.boundaryVerticesᶜ, M.angleDefect v = 0) (hχ : M.eulerChar = 1) (w : ℝ)
    {ℓ lam : {v // v ∈ M.boundaryVertices} → ℝ}
    {κ : {v // v ∈ M.boundaryVertices} → {v // v ∈ M.boundaryVertices} → ℝ}
    (hκ : ∀ p q, κ p q = κ q p)
    (h : fluxMatrix ℓ κ *ᵥ lam = fun p => Ddg.fluxRhs (w * (π - M.angleSum p.1))) :
    ∑ p, ℓ p * lam p = 2 * π * w := by
  rw [M.onCut_force_sum_const hM hman w hκ h, Finset.sum_eq_zero hflat, hχ]
  push_cast
  ring

end Turning

/-! ## The derivative of the weight across a segment -/

section Foot

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The foot of the height of `a` over the line `pq` is `p + t (q − p)` with
`t = ⟨a − p, q − p⟩ / |q − p|²`. -/
theorem inner_sub_foot (a p q : E) (hpq : p ≠ q) :
    inner ℝ (a - (p + (inner ℝ (a - p) (q - p) / ‖q - p‖ ^ 2) • (q - p))) (q - p) = 0 := by
  have hn : ‖q - p‖ ^ 2 ≠ 0 := by
    have : q - p ≠ 0 := sub_ne_zero.mpr (Ne.symm hpq)
    positivity
  have e : a - (p + (inner ℝ (a - p) (q - p) / ‖q - p‖ ^ 2) • (q - p)) =
      (a - p) - (inner ℝ (a - p) (q - p) / ‖q - p‖ ^ 2) • (q - p) := by abel
  rw [e, inner_sub_left, real_inner_smul_left, real_inner_self_eq_norm_sq,
    div_mul_cancel₀ _ hn, sub_self]

/-- The two cotangent numerators of the triangle `a p q` at `p` and `q` add up to `|pq|²`, so
the fraction `t` of the foot is `cot_p/(cot_p + cot_q)` (both cotangents are the numerators
divided by twice the area), as SPEC §3.10 writes it. -/
theorem inner_add_inner_eq_norm_sq (a p q : E) :
    inner ℝ (a - p) (q - p) + inner ℝ (a - q) (p - q) = ‖q - p‖ ^ 2 := by
  have e : p - q = -(q - p) := by abel
  rw [e, inner_neg_right, ← sub_eq_add_neg, ← inner_sub_left, ← real_inner_self_eq_norm_sq]
  congr 1
  abel

end Foot

end Cuts

end DdgProofs
