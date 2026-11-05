import DdgProofs.Cuts.Recovery

/-!
# Cut flow: the default stages on the cut

`crates/ddg-variational-cuts/SPEC.md` §3.8 – §3.10, D14 – D16: the length force and the rate are
computed on the cut, the rate is carried to the vertices by weighted means, and it is scaled by
the mean of `|∇φ|` over the piece of the cut.

The index-set definitions of this file (`wmean`, `localExtension`, `meanGradient`,
`rippleFactor`, `rippleRate`) are **this package's own**; the two linear systems are
`Cuts.fluxMatrix`, built from the oracle's `fluxDiagEntry`, `fluxOffEntry`, `fluxRhs`
(`DdgProofs/Cuts/Recovery.lean`). The last sections restate the results about the oracle's
scalar definitions `gradientNormLone`, `meanGradient`, `weightedMean`, `ringMean`,
`rateCoupling`, `rateRhs`, `segmentMeanWeight`, `turningAngle`, `lengthForce`, `angleFromCot`,
`footFraction`, `heightOver`, `normalDerivativeOfLinear`, `adaptiveStep`.

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

namespace DDG

namespace Cuts

open Matrix
open scoped DDG.Oracle
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

/-- **One face of small slope takes the harmonic mean of its piece with it** (SPEC D30): for a
face `σ₀` of the piece with `L_{σ₀} > 0`, `ḡ ≤ (ΣL / L_{σ₀}) · g_{σ₀}` (lengths `≥ 0`, slopes
`> 0`). So a slope that the formula gets wrong by a factor makes the mean wrong by that factor
times at most the number `ΣL / L_{σ₀}`; this is why a face with a moved crossing is left out.

[formula: cut-gradient-norm] -/
theorem meanGradient_le_of_face {ι : Type*} {s : Finset ι} {L g : ι → ℝ}
    (hL : ∀ σ ∈ s, 0 ≤ L σ) (hg : ∀ σ ∈ s, 0 < g σ) {σ₀ : ι} (h0 : σ₀ ∈ s) (hL0 : 0 < L σ₀) :
    meanGradient s L g ≤ (∑ σ ∈ s, L σ) / L σ₀ * g σ₀ := by
  have hsum : 0 ≤ ∑ σ ∈ s, L σ := Finset.sum_nonneg hL
  have hg0 := hg σ₀ h0
  have hterm : L σ₀ / g σ₀ ≤ ∑ σ ∈ s, L σ / g σ :=
    Finset.single_le_sum (f := fun σ => L σ / g σ)
      (fun σ hσ => div_nonneg (hL σ hσ) (hg σ hσ).le) h0
  have hpos : 0 < L σ₀ / g σ₀ := div_pos hL0 hg0
  have e : (∑ σ ∈ s, L σ) / L σ₀ * g σ₀ = (∑ σ ∈ s, L σ) / (L σ₀ / g σ₀) := by
    field_simp
  rw [e]
  unfold meanGradient
  exact div_le_div_of_nonneg_left hsum hpos hterm

/-- The mean slope of a piece as the flow uses it (SPEC §3.9, D23, D30): the harmonic mean over
the faces of the piece that count if that is positive, and `1` otherwise. `s` is the set of the
faces that count: any set, in the code the crossed faces none of whose crossings was moved and
whose slope is positive. -/
def pieceSlope {ι : Type*} (s : Finset ι) (L g : ι → ℝ) : ℝ :=
  if 0 < meanGradient s L g then meanGradient s L g else 1

/-- The scale `ḡ` of a piece is positive, whatever counts. -/
theorem pieceSlope_pos {ι : Type*} (s : Finset ι) (L g : ι → ℝ) : 0 < pieceSlope s L g := by
  unfold pieceSlope
  split_ifs with h
  · exact h
  · exact one_pos

/-- **A piece without any face that counts has `ḡ = 1`.**

[formula: cut-gradient-norm] -/
theorem pieceSlope_empty {ι : Type*} (L g : ι → ℝ) : pieceSlope ∅ L g = 1 := by
  unfold pieceSlope meanGradient
  simp

/-- Where faces count (lengths `≥ 0` with a positive sum, slopes `> 0`), `ḡ` is their harmonic
mean, to which `meanGradient_mem` and `meanGradient_normalizes` apply. -/
theorem pieceSlope_eq_meanGradient {ι : Type*} {s : Finset ι} {L g : ι → ℝ}
    (hL : ∀ σ ∈ s, 0 ≤ L σ) (hpos : 0 < ∑ σ ∈ s, L σ) (hg : ∀ σ ∈ s, 0 < g σ) :
    pieceSlope s L g = meanGradient s L g := by
  obtain ⟨σ, hσ, hlt⟩ : ∃ σ ∈ s, 0 < L σ := by
    by_contra h
    push Not at h
    have := Finset.sum_nonpos h
    linarith
  have hden : 0 < ∑ σ ∈ s, L σ / g σ :=
    Finset.sum_pos' (fun τ hτ => div_nonneg (hL τ hτ) (hg τ hτ).le)
      ⟨σ, hσ, div_pos hlt (hg σ hσ)⟩
  unfold pieceSlope
  rw [ite_eq_left]
  unfold meanGradient
  exact div_pos hpos hden

/-- **The faces that count after D30**: of the crossed faces `s` of a piece, those without a
moved crossing and with a positive slope. If every face has a moved crossing or the slope `0`,
none counts and `ḡ = 1`: with D30 a piece all of whose faces have a moved crossing is one more
way for this to happen.

[formula: cut-gradient-norm] -/
theorem pieceSlope_of_no_counted {ι : Type*} {s : Finset ι} (moved : ι → Prop)
    [DecidablePred moved] {L g : ι → ℝ} (h : ∀ σ ∈ s, moved σ ∨ ¬ 0 < g σ) :
    pieceSlope (s.filter fun σ => ¬ moved σ ∧ 0 < g σ) L g = 1 := by
  have e : s.filter (fun σ => ¬ moved σ ∧ 0 < g σ) = ∅ := by
    rw [Finset.filter_eq_empty_iff]
    intro σ hσ hc
    rcases h σ hσ with h1 | h1
    · exact hc.1 h1
    · exact h1 hc.2
  rw [e, pieceSlope_empty]

/-- On the faces that count the slopes are positive by the choice of the set, so the mean lies
between the smallest and the largest slope of the faces that count. -/
theorem pieceSlope_counted_mem {ι : Type*} {s : Finset ι} (moved : ι → Prop)
    [DecidablePred moved] {L g : ι → ℝ} (hL : ∀ σ ∈ s, 0 ≤ L σ)
    (hpos : 0 < ∑ σ ∈ s.filter (fun σ => ¬ moved σ ∧ 0 < g σ), L σ) {lo hi : ℝ} (hlo0 : 0 < lo)
    (hlo : ∀ σ ∈ s, ¬ moved σ → 0 < g σ → lo ≤ g σ)
    (hhi : ∀ σ ∈ s, ¬ moved σ → 0 < g σ → g σ ≤ hi) :
    lo ≤ pieceSlope (s.filter fun σ => ¬ moved σ ∧ 0 < g σ) L g ∧
      pieceSlope (s.filter fun σ => ¬ moved σ ∧ 0 < g σ) L g ≤ hi := by
  have hL' : ∀ σ ∈ s.filter (fun σ => ¬ moved σ ∧ 0 < g σ), 0 ≤ L σ :=
    fun σ hσ => hL σ (Finset.mem_filter.mp hσ).1
  have hg' : ∀ σ ∈ s.filter (fun σ => ¬ moved σ ∧ 0 < g σ), 0 < g σ :=
    fun σ hσ => (Finset.mem_filter.mp hσ).2.2
  rw [pieceSlope_eq_meanGradient hL' hpos hg']
  exact meanGradient_mem hL' hpos hlo0
    (fun σ hσ => hlo σ (Finset.mem_filter.mp hσ).1 (Finset.mem_filter.mp hσ).2.1
      (Finset.mem_filter.mp hσ).2.2)
    (fun σ hσ => hhi σ (Finset.mem_filter.mp hσ).1 (Finset.mem_filter.mp hσ).2.1
      (Finset.mem_filter.mp hσ).2.2)

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
theorem _root_.DDG.TriMesh.sum_turningAngle (hM : M.Valid) (hman : M.Manifold) :
    ∑ v ∈ M.boundaryVertices, (π - M.angleSum v) =
      2 * π * (M.eulerChar : ℝ) - ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v := by
  have h := M.gauss_bonnet hM hman
  simp only [TriMesh.geodesicCurvature] at h
  linarith

/-- **A closed cut on a flat surface**: if the region is a flat disk (no angle defect at its
interior vertices, Euler characteristic `1`), its boundary turns by `2π` in total. -/
theorem _root_.DDG.TriMesh.sum_turningAngle_flat_disk (hM : M.Valid) (hman : M.Manifold)
    (hflat : ∀ v ∈ M.boundaryVerticesᶜ, M.angleDefect v = 0) (hχ : M.eulerChar = 1) :
    ∑ v ∈ M.boundaryVertices, (π - M.angleSum v) = 2 * π := by
  rw [M.sum_turningAngle hM hman, Finset.sum_eq_zero hflat, hχ]
  push_cast
  ring

/-- **The total length force of a closed cut with a constant weight**: with
`Θ_p = w (π − β_p)` at the boundary vertices of the cut mesh of the region,
`Σ_p ℓ_p λ⁰_p = w (2πχ − Σ_interior Ω)`.

[formula: cut-length-step-on-cut] -/
theorem _root_.DDG.TriMesh.onCut_force_sum_const (hM : M.Valid) (hman : M.Manifold)
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
theorem _root_.DDG.TriMesh.onCut_force_sum_flat_disk (hM : M.Valid) (hman : M.Manifold)
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

/-! ## The oracle's definitions

`lean/oracle/DdgOracle/Formulas/Cuts.lean` (revision "default stages") defines the scalar
pieces of these stages. The theorems above are restated about them here. -/

/-- The oracle's `gradientNormLone` is `gradNorm`. -/
theorem gradientNormLone_eq (φa len area : ℝ) :
    Ddg.gradientNormLone φa len area = gradNorm φa len area := by
  simp only [Ddg.gradientNormLone, gradNorm, Oracle.abs_eq, Oracle.two_eq]

/-- **`|∇φ|` in a crossed face, on the oracle's definition**: for the gradient
`g = α e₁ + β e₂` of an affine function with the value `φ_a` at `a` that vanishes at
`p = a + e₁` and `q = a + e₂`, `‖g‖ = gradientNormLone φ_a |pq| area(a, p, q)`.

[formula: cut-gradient-norm] -/
theorem norm_grad_eq_gradientNormLone {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] (e₁ e₂ : E) (α β φa : ℝ) (h : IsTriangle ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖)
    (h1 : inner ℝ (α • e₁ + β • e₂) e₁ = -φa) (h2 : inner ℝ (α • e₁ + β • e₂) e₂ = -φa) :
    ‖α • e₁ + β • e₂‖ = Ddg.gradientNormLone φa ‖e₂ - e₁‖ (triArea ‖e₂ - e₁‖ ‖e₂‖ ‖e₁‖) := by
  rw [gradientNormLone_eq]
  exact norm_grad_eq_gradNorm e₁ e₂ α β φa h h1 h2

/-- **A zero-slope face** (SPEC §3.9, D23a): if the function vanishes at the lone vertex,
`|∇φ| = 0` in the face. Such a face is left out of the mean `ḡ` (its term `|pq|/|∇φ|` is not
finite); `meanGradient_mem` then applies to the faces that count, where `|∇φ| > 0`. -/
theorem gradientNormLone_zero (len area : ℝ) : Ddg.gradientNormLone 0 len area = 0 := by
  rw [gradientNormLone_eq]
  simp [gradNorm]

/-- The oracle's `meanGradient` of the two sums is `meanGradient`. -/
theorem oracle_meanGradient_eq {ι : Type*} (s : Finset ι) (L g : ι → ℝ) :
    Ddg.meanGradient (∑ σ ∈ s, L σ) (∑ σ ∈ s, L σ / g σ) = meanGradient s L g := rfl

/-- **The oracle's `ḡ` lies between the smallest and the largest `|∇φ|` of the piece.**

[formula: cut-gradient-norm] -/
theorem oracle_meanGradient_mem {ι : Type*} {s : Finset ι} {L g : ι → ℝ}
    (hL : ∀ σ ∈ s, 0 ≤ L σ) (hpos : 0 < ∑ σ ∈ s, L σ) {lo hi : ℝ} (hlo0 : 0 < lo)
    (hlo : ∀ σ ∈ s, lo ≤ g σ) (hhi : ∀ σ ∈ s, g σ ≤ hi) :
    lo ≤ Ddg.meanGradient (∑ σ ∈ s, L σ) (∑ σ ∈ s, L σ / g σ) ∧
      Ddg.meanGradient (∑ σ ∈ s, L σ) (∑ σ ∈ s, L σ / g σ) ≤ hi := by
  rw [oracle_meanGradient_eq]
  exact meanGradient_mem hL hpos hlo0 hlo hhi

/-- The oracle's `weightedMean` at `ℝ`: the two lists are zipped, and the sum of the products
is divided by the sum of the weights that are used. -/
theorem weightedMean_eq (ws vs : List ℝ) :
    Ddg.weightedMean ws vs =
      ((ws.zip vs).map fun y => y.1 * y.2).sum / ((ws.zip vs).map fun y => y.1).sum := by
  show Ddg.sumList ((ws.zip vs).map fun y => y.1 * y.2) /
      Ddg.sumList ((ws.zip vs).map fun y => y.1) = _
  rw [sumList_eq_sum, sumList_eq_sum]

/-- **The oracle's `weightedMean` is a convex combination** of the values it uses: for
nonnegative weights with a positive sum it lies between every lower and upper bound of the values
(of the pairs of the two lists, whatever their lengths).

[formula: cut-extension-local] -/
theorem weightedMean_mem {ws vs : List ℝ} (hw : ∀ y ∈ ws.zip vs, 0 ≤ y.1)
    (hpos : 0 < ((ws.zip vs).map fun y => y.1).sum) {lo hi : ℝ}
    (hlo : ∀ y ∈ ws.zip vs, lo ≤ y.2) (hhi : ∀ y ∈ ws.zip vs, y.2 ≤ hi) :
    lo ≤ Ddg.weightedMean ws vs ∧ Ddg.weightedMean ws vs ≤ hi := by
  rw [weightedMean_eq]
  generalize ws.zip vs = ps at *
  have key : lo * (ps.map fun y => y.1).sum ≤ (ps.map fun y => y.1 * y.2).sum ∧
      (ps.map fun y => y.1 * y.2).sum ≤ hi * (ps.map fun y => y.1).sum := by
    clear hpos
    induction ps with
    | nil => simp
    | cons y t ih =>
      have h0 := hw y (List.mem_cons.mpr (Or.inl rfl))
      have h1 := hlo y (List.mem_cons.mpr (Or.inl rfl))
      have h2 := hhi y (List.mem_cons.mpr (Or.inl rfl))
      obtain ⟨ih1, ih2⟩ := ih (fun z hz => hw z (List.mem_cons.mpr (Or.inr hz)))
        (fun z hz => hlo z (List.mem_cons.mpr (Or.inr hz)))
        (fun z hz => hhi z (List.mem_cons.mpr (Or.inr hz)))
      simp only [List.map_cons, List.sum_cons]
      constructor <;> nlinarith [mul_le_mul_of_nonneg_left h1 h0, mul_le_mul_of_nonneg_left h2 h0]
  rw [le_div_iff₀ hpos, div_le_iff₀ hpos]
  exact key

/-- **The oracle's `weightedMean` of two lists is the weighted mean** over the (distinct) indices
the lists are built from. -/
theorem oracle_weightedMean_eq {ι : Type*} [DecidableEq ι] {l : List ι} (hl : l.Nodup)
    (c x : ι → ℝ) : Ddg.weightedMean (l.map c) (l.map x) = wmean l.toFinset c x := by
  have h1 : ∀ l : List ι,
      (((l.map c).zip (l.map x)).map fun y => y.1 * y.2).sum =
        (l.map fun j => c j * x j).sum := by
    intro l
    induction l with
    | nil => simp
    | cons a t ih => simp [ih]
  have h2 : ∀ l : List ι, (((l.map c).zip (l.map x)).map fun y => y.1).sum = (l.map c).sum := by
    intro l
    induction l with
    | nil => simp
    | cons a t ih => simp [ih]
  rw [weightedMean_eq, h1, h2]
  unfold wmean
  rw [← List.sum_toFinset (fun j => c j * x j) hl, ← List.sum_toFinset c hl]

/-- **The first line of the local extension, on the oracle's definition**: `weightedMean` of
the weights `ℓ_p β_ip ≥ 0` and the rates `r_p` over the cut points at `i` lies between the
smallest and the largest of these rates.

[formula: cut-extension-local] -/
theorem oracle_weightedMean_mem {ι : Type*} [DecidableEq ι] {l : List ι} (hl : l.Nodup)
    {c x : ι → ℝ} (hc : ∀ j ∈ l, 0 ≤ c j) (hpos : 0 < ∑ j ∈ l.toFinset, c j) {lo hi : ℝ}
    (hlo : ∀ j ∈ l, lo ≤ x j) (hhi : ∀ j ∈ l, x j ≤ hi) :
    lo ≤ Ddg.weightedMean (l.map c) (l.map x) ∧ Ddg.weightedMean (l.map c) (l.map x) ≤ hi := by
  rw [oracle_weightedMean_eq hl]
  exact wmean_mem (fun j hj => hc j (List.mem_toFinset.mp hj)) hpos
    (fun j hj _ => hlo j (List.mem_toFinset.mp hj)) (fun j hj _ => hhi j (List.mem_toFinset.mp hj))

/-- **The second line of the local extension, on the oracle's definition**: `ringMean` of the
pairs `(σ_j, ℓ_ij)` of the neighbours one ring closer (positive edge lengths, at least one
neighbour) lies between the smallest and the largest `σ_j`.

[formula: cut-extension-local] -/
theorem oracle_ringMean_mem {ι : Type*} [DecidableEq ι] {l : List ι} (hl : l.Nodup)
    (hne : l ≠ []) {σ len : ι → ℝ} (hlen : ∀ j ∈ l, 0 < len j) {lo hi : ℝ}
    (hlo : ∀ j ∈ l, lo ≤ σ j) (hhi : ∀ j ∈ l, σ j ≤ hi) :
    lo ≤ Ddg.ringMean (l.map fun j => (σ j, len j)) ∧
      Ddg.ringMean (l.map fun j => (σ j, len j)) ≤ hi := by
  have e : Ddg.ringMean (l.map fun j => (σ j, len j)) =
      (∑ j ∈ l.toFinset, σ j / len j) / ∑ j ∈ l.toFinset, 1 / len j := by
    unfold Ddg.ringMean
    simp only [List.map_map, sumList_eq_sum, Oracle.one_eq, Function.comp_def]
    rw [← List.sum_toFinset (fun j => σ j / len j) hl, ← List.sum_toFinset (fun j => 1 / len j) hl]
  rw [e]
  have hN : l.toFinset.Nonempty := by
    obtain ⟨a, ha⟩ := List.exists_mem_of_ne_nil l hne
    exact ⟨a, List.mem_toFinset.mpr ha⟩
  exact outwardExtension_mem hN (fun j hj => hlen j (List.mem_toFinset.mp hj))
    (fun j hj => hlo j (List.mem_toFinset.mp hj)) (fun j hj => hhi j (List.mem_toFinset.mp hj))

/-- The oracle's `rateCoupling`, `rateRhs`, `segmentMeanWeight` at `ℝ`. -/
theorem rateCoupling_eq (τ wbar len : ℝ) : Ddg.rateCoupling τ wbar len = τ * (wbar / len) := by
  unfold Ddg.rateCoupling
  ring

theorem rateRhs_eq (ℓ s lam : ℝ) : Ddg.rateRhs ℓ s lam = ℓ * (s + lam) := rfl

theorem segmentMeanWeight_eq (wp wq : ℝ) : Ddg.segmentMeanWeight wp wq = (wp + wq) / 2 := by
  simp only [Ddg.segmentMeanWeight, Oracle.two_eq]

theorem segmentMeanWeight_nonneg {wp wq : ℝ} (hp : 0 ≤ wp) (hq : 0 ≤ wq) :
    0 ≤ Ddg.segmentMeanWeight wp wq := by
  rw [segmentMeanWeight_eq]
  positivity

/-- **The rate system on the oracle's entries** (SPEC §3.10): couplings
`rateCoupling τ w̄_σ |σ| = τ w̄_σ/|σ|` (symmetric, `w̄ ≥ 0`, `|σ| ≥ 0`, `τ ≥ 0`), right-hand side
`rateRhs ℓ_p s_p λ⁰_p = ℓ_p (s_p + λ⁰_p)`, positive `ℓ`. Its solution `r` satisfies
`Σ ℓ r = Σ ℓ (s + λ⁰)`, lies between every lower and upper bound of `s + λ⁰`, and vanishes if
`s + λ⁰ = 0`, whatever `τ`.

[formula: cut-length-step-on-cut] -/
theorem onCut_rate_oracle {C : Type*} [Fintype C] [DecidableEq C] {ℓ s lam r : C → ℝ}
    {wbar len : C → C → ℝ} {τ : ℝ} (hℓ : ∀ p, 0 < ℓ p) (hw : ∀ p q, wbar p q = wbar q p)
    (hw0 : ∀ p q, 0 ≤ wbar p q) (hlen : ∀ p q, len p q = len q p) (hlen0 : ∀ p q, 0 ≤ len p q)
    (hτ : 0 ≤ τ)
    (h : fluxMatrix ℓ (fun p q => Ddg.rateCoupling τ (wbar p q) (len p q)) *ᵥ r =
      fun p => Ddg.rateRhs (ℓ p) (s p) (lam p)) :
    (∑ p, ℓ p * r p = ∑ p, ℓ p * (s p + lam p)) ∧
      (∀ lo hi : ℝ, (∀ p, lo ≤ s p + lam p) → (∀ p, s p + lam p ≤ hi) →
        ∀ p, lo ≤ r p ∧ r p ≤ hi) ∧
      ((∀ p, s p + lam p = 0) → ∀ p, r p = 0) := by
  have h' : fluxMatrix ℓ (fun p q => τ * (wbar p q / len p q)) *ᵥ r =
      fun p => Ddg.fluxRhs (ℓ p * (s p + lam p)) := by
    simpa only [rateCoupling_eq, rateRhs_eq, fluxRhs_eq] using h
  have hm : ∀ p q, wbar p q / len p q = wbar q p / len q p := fun p q => by
    rw [hw p q, hlen p q]
  have hm0 : ∀ p q, 0 ≤ wbar p q / len p q := fun p q => div_nonneg (hw0 p q) (hlen0 p q)
  exact ⟨onCut_rate_sum hm h', fun lo hi hlo hhi p => onCut_rate_bounds hℓ hm0 hτ h' hlo hhi p,
    fun h0 p => onCut_rest hℓ hm0 hτ h' h0 p⟩

/-- The oracle's `turningAngle` and `lengthForce` at `ℝ`. -/
theorem turningAngle_eq (γ β : ℝ) : Ddg.turningAngle γ β = γ - β := rfl

theorem lengthForce_eq (w γ β : ℝ) (segs : List (ℝ × ℝ)) :
    Ddg.lengthForce w γ β segs = w * (γ - β) + (segs.map fun x => x.1 * x.2 / 2).sum := by
  simp only [Ddg.lengthForce, turningAngle_eq, sumList_eq_sum, Oracle.two_eq]

/-- With a constant weight (`∂w/∂n = 0` on every segment) the length force is `w (γ − β)`. -/
theorem lengthForce_of_const (w γ β : ℝ) {segs : List (ℝ × ℝ)} (h : ∀ x ∈ segs, x.2 = 0) :
    Ddg.lengthForce w γ β segs = w * (γ - β) := by
  rw [lengthForce_eq]
  have : (segs.map fun x => x.1 * x.2 / 2).sum = 0 := by
    apply List.sum_eq_zero
    intro y hy
    obtain ⟨x, hx, rfl⟩ := List.mem_map.mp hy
    rw [h x hx]
    ring
  rw [this, add_zero]

/-- **The total length force of a closed cut, on the oracle's `lengthForce`**: with
`Θ_p = lengthForce w π β_p segs_p` at the boundary vertices of the cut mesh of a region, `β_p`
the sum of the corner angles there and a constant weight `w`,
`Σ_p ℓ_p λ⁰_p = w (2πχ − Σ_interior Ω)`; for a flat disk this is `2πw`.

[formula: cut-length-step-on-cut] -/
theorem _root_.DDG.TriMesh.onCut_lengthForce_sum {V F : Type*} [DecidableEq V] [Fintype V]
    [Fintype F] (M : TriMesh V F) (hM : M.Valid) (hman : M.Manifold) (w : ℝ)
    {ℓ lam : {v // v ∈ M.boundaryVertices} → ℝ}
    {κ : {v // v ∈ M.boundaryVertices} → {v // v ∈ M.boundaryVertices} → ℝ}
    (hκ : ∀ p q, κ p q = κ q p) (segs : {v // v ∈ M.boundaryVertices} → List (ℝ × ℝ))
    (hsegs : ∀ p, ∀ x ∈ segs p, x.2 = 0)
    (h : fluxMatrix ℓ κ *ᵥ lam =
      fun p => Ddg.fluxRhs (Ddg.lengthForce w π (M.angleSum p.1) (segs p))) :
    ∑ p, ℓ p * lam p =
      w * (2 * π * (M.eulerChar : ℝ) - ∑ v ∈ M.boundaryVerticesᶜ, M.angleDefect v) := by
  refine M.onCut_force_sum_const hM hman w hκ ?_
  rw [h]
  funext p
  rw [lengthForce_of_const w π _ (hsegs p)]

/-- The oracle's `angleFromCot` recovers an angle of `(0, π)` from its cotangent. -/
theorem angleFromCot_cot {θ : ℝ} (h0 : 0 < θ) (hπ : θ < π) :
    Ddg.angleFromCot (Real.cot θ) = θ := by
  simp only [Ddg.angleFromCot, Oracle.pi_eq, Oracle.two_eq, Oracle.atan_eq]
  have e : Real.cot θ = Real.tan (π / 2 - θ) := by
    rw [Real.tan_pi_div_two_sub, Real.tan_eq_sin_div_cos, Real.cot_eq_cos_div_sin, inv_div]
  rw [e, Real.arctan_tan (by linarith) (by linarith)]
  ring

/-- **The interior angle from the cotangents of the cut mesh**: `angleFromCot` of the cotangent
of a corner of a triangle is its corner angle, so the sum of `angleFromCot` over the triangles of
a region at a boundary node is the interior angle `β` of SPEC §3.10.

[formula: cut-length-step-on-cut] -/
theorem angleFromCot_cornerCot {a b c : ℝ} (h : IsTriangle a b c) :
    Ddg.angleFromCot (cornerCot a b c) = cornerAngle a b c := by
  rw [cornerCot_eq_cot h]
  exact angleFromCot_cot (cornerAngle_pos h) (cornerAngle_lt_pi h)

/-- The oracle's `footFraction` in a triangle `z p q` with side lengths `lzp`, `lzq`, `lpq`: from
the cotangents at `p` and `q` it is `(lzp² + lpq² − lzq²)/(2 lpq²)`. -/
theorem footFraction_eq {lzp lzq lpq : ℝ} (h : IsTriangle lzq lzp lpq) :
    Ddg.footFraction (cornerCot lzq lzp lpq) (cornerCot lzp lpq lzq) =
      (lzp ^ 2 + lpq ^ 2 - lzq ^ 2) / (2 * lpq ^ 2) := by
  have hT := (triArea_pos h).ne'
  have hpq := h.pos_c.ne'
  unfold Ddg.footFraction cornerCot
  rw [triArea_rotate lzq lzp lpq]
  have hsum : (lzp ^ 2 + lpq ^ 2 - lzq ^ 2) / (4 * triArea lzq lzp lpq) +
      (lpq ^ 2 + lzq ^ 2 - lzp ^ 2) / (4 * triArea lzq lzp lpq) =
      2 * lpq ^ 2 / (4 * triArea lzq lzp lpq) := by ring
  rw [hsum, div_div_div_cancel_right₀ (by positivity)]

/-- **The foot of the height, from the cotangents** (SPEC §3.10): for points `z, p, q` of a real
inner product space forming a nondegenerate triangle, the point at the oracle's `footFraction`
of the cotangents at `p` and `q` along `pq` is the foot of the height from `z`: the vector from
it to `z` is orthogonal to `pq`.

[formula: cut-length-step-on-cut] -/
theorem footFraction_orthogonal {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (z p q : E) (h : IsTriangle ‖z - q‖ ‖z - p‖ ‖q - p‖) :
    inner ℝ (z - (p + Ddg.footFraction (cornerCot ‖z - q‖ ‖z - p‖ ‖q - p‖)
      (cornerCot ‖z - p‖ ‖q - p‖ ‖z - q‖) • (q - p))) (q - p) = 0 := by
  have hpq : p ≠ q := by
    intro hpq
    have := h.pos_c
    rw [hpq, sub_self, norm_zero] at this
    exact lt_irrefl _ this
  have e : Ddg.footFraction (cornerCot ‖z - q‖ ‖z - p‖ ‖q - p‖)
      (cornerCot ‖z - p‖ ‖q - p‖ ‖z - q‖) = inner ℝ (z - p) (q - p) / ‖q - p‖ ^ 2 := by
    rw [footFraction_eq h]
    have h1 : z - q = (z - p) - (q - p) := by abel
    have h2 : ‖z - q‖ ^ 2 = ‖z - p‖ ^ 2 - 2 * inner ℝ (z - p) (q - p) + ‖q - p‖ ^ 2 := by
      rw [h1, norm_sub_sq_real]
    have hn : ‖q - p‖ ^ 2 ≠ 0 := by
      have := h.pos_c
      positivity
    rw [h2]
    field_simp
    ring
  rw [e]
  exact inner_sub_foot z p q hpq

/-- The oracle's `normalDerivativeOfLinear` for an affine function `w(x) = w₀ + ⟨g, x⟩`: the
increment of `w` from the third vertex `z` to the point `f = (1 − λ) p + λ q` of the segment,
divided by `h`. With `λ` the foot fraction and `h = heightOver area |pq|` the height, this is the
derivative of `w` in the direction from `z` across the segment. -/
theorem normalDerivativeOfLinear_eq {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (g z p q : E) (w0 lam h : ℝ) :
    Ddg.normalDerivativeOfLinear (w0 + inner ℝ g p) (w0 + inner ℝ g q) (w0 + inner ℝ g z) lam h =
      inner ℝ g (((1 - lam) • p + lam • q) - z) / h := by
  simp only [Ddg.normalDerivativeOfLinear, Oracle.one_eq, inner_sub_right, inner_add_right,
    real_inner_smul_right]
  congr 1
  ring

/-- The oracle's `heightOver`: `2·area/|pq|`, the height of the third vertex over the segment
(`height · |pq| = 2 · area`). -/
theorem heightOver_mul (area : ℝ) {len : ℝ} (hlen : len ≠ 0) :
    Ddg.heightOver area len * len = 2 * area := by
  simp only [Ddg.heightOver, Oracle.two_eq]
  field_simp

/-! ## The adaptive step (SPEC §4)

`adaptiveStep τ_max time_scale ω max_motion h ρ̂` with the statistics `ω`, `h`, `ρ̂` of the cut:
the largest local rate, the shortest crossed edge and the largest `|ρ_c|` for the extremes;
the quantiles `ω`, `m̂` and `h = 1` for the robust statistics (the oracle's `robustStep`,
`DdgProofs/Cuts/Quantile.lean`, which also has the skipping of terms and the local limit). -/

/-- The oracle's `adaptiveStep` at `ℝ`. -/
theorem adaptiveStep_eq (τmax ts ω mm h ρ : ℝ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ =
      if 0 < ρ then min (if 0 < ω then min τmax (ts / ω) else τmax) (mm * h / ρ)
      else (if 0 < ω then min τmax (ts / ω) else τmax) := by
  simp only [Ddg.adaptiveStep, Oracle.zero_eq]

/-- **The adaptive step is at most `τ_max`.**

[formula: cut-step-cap] -/
theorem adaptiveStep_le_max (τmax ts ω mm h ρ : ℝ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ ≤ τmax := by
  rw [adaptiveStep_eq]
  split_ifs
  · exact (min_le_left _ _).trans (min_le_left _ _)
  · exact min_le_left _ _
  · exact min_le_left _ _
  · exact le_rfl

/-- **The adaptive step resolves the roughening rate**: `τ · ω ≤ time_scale` for `ω > 0`.

[formula: cut-step-cap] -/
theorem adaptiveStep_mul_omega_le (τmax ts mm h ρ : ℝ) {ω : ℝ} (hω : 0 < ω) :
    Ddg.adaptiveStep τmax ts ω mm h ρ * ω ≤ ts := by
  have h1 : Ddg.adaptiveStep τmax ts ω mm h ρ ≤ ts / ω := by
    rw [adaptiveStep_eq]
    split_ifs
    · exact (min_le_left _ _).trans (min_le_right _ _)
    · exact min_le_right _ _
  calc Ddg.adaptiveStep τmax ts ω mm h ρ * ω ≤ ts / ω * ω :=
        mul_le_mul_of_nonneg_right h1 hω.le
    _ = ts := by field_simp

/-- **The adaptive step limits the motion**: `τ · ρ̂ ≤ max_motion · h` for a speed statistic
`ρ̂ > 0`.

[formula: cut-step-cap] -/
theorem adaptiveStep_mul_rho_le (τmax ts ω mm h : ℝ) {ρ : ℝ} (hρ : 0 < ρ) :
    Ddg.adaptiveStep τmax ts ω mm h ρ * ρ ≤ mm * h := by
  have h1 : Ddg.adaptiveStep τmax ts ω mm h ρ ≤ mm * h / ρ := by
    rw [adaptiveStep_eq, ite_eq_left hρ]
    exact min_le_right _ _
  calc Ddg.adaptiveStep τmax ts ω mm h ρ * ρ ≤ mm * h / ρ * ρ :=
        mul_le_mul_of_nonneg_right h1 hρ.le
    _ = mm * h := by field_simp

end Cuts

end DDG
