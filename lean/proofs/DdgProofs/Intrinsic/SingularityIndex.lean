import DdgProofs.Orientation
import DdgProofs.Intrinsic.AppsAmrVector

/-!
# The index of a face and the Poincaré–Hopf sum (`app-vector-singularity-index`)

`vector::singularities` of `crates/ddg-intrinsic-apps/src/vector.rs`, on a closed oriented
valid mesh, in exact arithmetic. The model, halfedge `h = (f, k)` from corner `k` to `k + 1`:

* `φ h`: the direction of `h` at its tail in the tangent frame of the tail (`tail_angle`);
  the frames are **signposts**: at every corner the direction rotates from the outgoing
  halfedge `(f, k + 1)` to the outgoing halfedge `twin (f, k)` by the scaled corner angle
  `s θ` (up to full turns, `SignpostFrames`), `s_v = 2π/Θ_v` (`tangent_scale`, interior);
* `ρ h = φ(twin h) + π − φ h` (`transport_angle`, VHM Eq. 10; on a closed mesh
  `head_back_angle` is the twin's tail angle), `r_h = e^{i n ρ_h}` (`transport`);
* `Ω_f = Σ_k (s_{c_k} − 1) θ_k` (`face_holonomy`);
* `ω_h = arg(u_j / (r_h u_i)) ∈ (−π, π]` for `h = i → j` (`(uj / (r * ui)).arg()`; Mathlib's
  `Complex.arg` has the same range as `f64::atan2`);
* `index f = (Σ_k ω_(f,k) + n Ω_f) / 2π`.

Proved, for `ω` evaluated by `arg` on **both** halfedges (the code before the fix):

* `faceIndex_isInt`: every `index f` is an integer (so the code's `round` returns it);
* `omega_twin`: `ω_{twin h} = −ω_h`, **except** `ω_{twin h} = π` when `ω_h = π`;
* `sum_faceIndex`, `sum_round_faceIndex`: `Σ_f index f = n χ + #{h : ω_h = π} / 2`: every edge
  at which the transported values are exactly opposite (`u_j = −c·r_ij u_i`, `c > 0`) adds `+1`;
* `sum_faceIndex_of_ne_pi`: if no `ω_h` equals `π`, `Σ_f index f = n χ`.

The code now (`singularities`) shifts `ω` on the halfedge with the higher index by `∓2π` when
the pair sum lies outside `[−π, π]` (`omegaCode`); in exact arithmetic that is `omegaFixed`:
`−ω_{twin h}` on the higher halfedge, at a tie `+π` on the lower and `−π` on the higher
(`omegaCode_eq_omegaFixed`). For that rule:

* `omegaFixed_add_twin`: the two halfedges of an edge have opposite angles, exactly;
* `faceIndexFixed_isInt`: every face index is still an integer (`omegaFixed ≡ ω` mod `2π`);
* `sum_faceIndexFixed`, `sum_round_faceIndexFixed`: **`index_sum = n χ` unconditionally** (discrete
  Poincaré–Hopf, GODF Eq. 19, App. B), on a closed valid oriented manifold mesh in which every
  vertex is a corner of some face (so `Θ_v > 0` and `s_v = 2π/Θ_v` is what the code computes;
  an isolated vertex would have `Θ_v = 0`).

`crates/ddg-intrinsic-apps/tests/vector_pi_tie.rs` is the regression test in floating point
(before the fix 25 of 240 single ties on icosphere(1) gave `index_sum = 3 ≠ 2`).
-/

noncomputable section

namespace DDG

open Real Complex

/-- Three nonzero complex numbers whose product is `e^{−iθ}`: their arguments and `θ` add up to
a multiple of `2π`. -/
theorem arg_add_arg_add_arg_eq {z₀ z₁ z₂ : ℂ} {θ : ℝ}
    (hz : z₀ * z₁ * z₂ * Complex.exp (θ * I) = 1) :
    ∃ m : ℤ, arg z₀ + arg z₁ + arg z₂ + θ = 2 * π * m := by
  have e : (‖z₀‖ * ‖z₁‖ * ‖z₂‖ : ℂ) * Complex.exp (((arg z₀ + arg z₁ + arg z₂ + θ : ℝ) : ℂ) * I) = 1 := by
    rw [← hz]
    conv_rhs => rw [← norm_mul_exp_arg_mul_I z₀, ← norm_mul_exp_arg_mul_I z₁,
      ← norm_mul_exp_arg_mul_I z₂]
    push_cast
    simp only [add_mul, Complex.exp_add]
    ring
  have hn : ‖z₀‖ * ‖z₁‖ * ‖z₂‖ = 1 := by
    have := congrArg (‖·‖) e
    simp only [norm_mul, Complex.norm_exp_ofReal_mul_I, mul_one, norm_one, Complex.norm_real,
      Real.norm_eq_abs, abs_norm] at this
    exact this
  have hn' : (‖z₀‖ * ‖z₁‖ * ‖z₂‖ : ℂ) = 1 := by exact_mod_cast hn
  rw [hn', one_mul, Complex.exp_eq_one_iff] at e
  obtain ⟨m, hm⟩ := e
  refine ⟨m, ?_⟩
  have := congrArg Complex.im hm
  simp at this
  linarith

namespace TriMesh

variable {V F : Type*} (M : TriMesh V F)

/-- The corner angle of `f` at corner `k` (`corner_angle(3f + k)`). -/
def cornerAngleAt (f : F) : Fin 3 → ℝ := ![M.angleA f, M.angleB f, M.angleC f]

variable [DecidableEq V] [Fintype V] [Fintype F]

/-- The scale from corner angles to tangent angles at an interior vertex, `s_v = 2π/Θ_v`
(`tangent_scale`). -/
def tangentScale (v : V) : ℝ := 2 * π / M.angleSum v

/-- The curvature the frames put into `f`, `Ω_f = Σ_k (s_{c_k} − 1) θ_k` (`face_holonomy`). -/
def faceHolonomy (f : F) : ℝ :=
  ∑ k, (M.tangentScale (M.corner f k) - 1) * M.cornerAngleAt f k

/-- Signpost frames: at corner `k + 1` of `f` the direction of the outgoing halfedge
`twin (f, k)` is that of the outgoing halfedge `(f, k + 1)` plus the scaled corner angle, up to
full turns. -/
def SignpostFrames (φ : F × Fin 3 → ℝ) : Prop :=
  ∀ f k, ∃ m : ℤ, φ (M.twin (f, k)) =
    φ (f, k + 1) + M.tangentScale (M.corner f (k + 1)) * M.cornerAngleAt f (k + 1) + 2 * π * m

/-- The transport angle `ρ_h = φ(twin h) + π − φ(h)` (`transport_angle`). -/
def transportAngle (φ : F × Fin 3 → ℝ) (h : F × Fin 3) : ℝ := φ (M.twin h) + π - φ h

/-- `ω_h = arg(u_j / (r_h u_i))` for `h = i → j`, `r_h = e^{i n ρ_h}`. -/
def omega (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ) (h : F × Fin 3) : ℝ :=
  arg (u (M.halfedge h).2 /
    (Complex.exp (((n * M.transportAngle φ h : ℝ) : ℂ) * I) * u (M.halfedge h).1))

/-- The index of `f` in units of `1/n`, `(Σ_k ω_(f,k) + n Ω_f)/2π` (GODF Eq. 19). -/
def faceIndex (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ) (f : F) : ℝ :=
  (∑ k, M.omega n φ u (f, k) + n * M.faceHolonomy f) / (2 * π)

omit [DecidableEq V] [Fintype V] [Fintype F] in
theorem transportAngle_add_twin (hO : M.Oriented) (φ : F × Fin 3 → ℝ) (h : F × Fin 3) :
    M.transportAngle φ h + M.transportAngle φ (M.twin h) = 2 * π := by
  simp only [transportAngle, M.twin_twin hO]
  ring

omit [DecidableEq V] [Fintype V] [Fintype F] in
/-- **The twin of a halfedge**: `ω_{twin h} = −ω_h`, except at `ω_h = π`, where both are `π`. -/
theorem omega_twin (hO : M.Oriented) (n : ℕ) (φ : F × Fin 3 → ℝ) {u : V → ℂ}
    (hu : ∀ v, u v ≠ 0) {h : F × Fin 3} (ht : M.HasTwin h) :
    M.omega n φ u (M.twin h) = if M.omega n φ u h = π then π else -M.omega n φ u h := by
  have hr : Complex.exp (((n * M.transportAngle φ (M.twin h) : ℝ) : ℂ) * I) =
      (Complex.exp (((n * M.transportAngle φ h : ℝ) : ℂ) * I))⁻¹ := by
    rw [← Complex.exp_neg]
    have e : M.transportAngle φ (M.twin h) = 2 * π - M.transportAngle φ h := by
      linarith [M.transportAngle_add_twin hO φ h]
    rw [e]
    have : (((n * (2 * π - M.transportAngle φ h) : ℝ) : ℂ) * I) =
        -(((n * M.transportAngle φ h : ℝ) : ℂ) * I) + (n : ℤ) * (2 * π * I) := by
      push_cast
      ring
    rw [this, Complex.exp_add, Complex.exp_int_mul_two_pi_mul_I, mul_one]
  unfold omega
  rw [M.halfedge_twin ht, Prod.fst_swap, Prod.snd_swap, hr, ← Complex.arg_inv]
  congr 1
  have h1 := hu (M.halfedge h).1
  have h2 := hu (M.halfedge h).2
  have h3 := Complex.exp_ne_zero (((n * M.transportAngle φ h : ℝ) : ℂ) * I)
  field_simp

omit [Fintype V] in
/-- The transport angles around a face add up to its holonomy, up to full turns
(`transport_around_face` in the halfedge model). -/
theorem sum_transportAngle (hM : M.Valid) {φ : F × Fin 3 → ℝ} (hφ : M.SignpostFrames φ) (f : F) :
    ∃ m : ℤ, ∑ k, M.transportAngle φ (f, k) = M.faceHolonomy f + 2 * π * m := by
  obtain ⟨m0, h0⟩ := hφ f 0
  obtain ⟨m1, h1⟩ := hφ f 1
  obtain ⟨m2, h2⟩ := hφ f 2
  have hsum := cornerAngle_add_cornerAngle_add_cornerAngle (hM.tri f)
  refine ⟨m0 + m1 + m2 + 2, ?_⟩
  simp only [Fin.sum_univ_three, transportAngle, faceHolonomy, cornerAngleAt, corner] at h0 h1 h2 ⊢
  simp only [Fin.isValue, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons, show (0 : Fin 3) + 1 = 1 from rfl,
    show (1 : Fin 3) + 1 = 2 from rfl, show (2 : Fin 3) + 1 = 0 from rfl] at h0 h1 h2 ⊢
  rw [h0, h1, h2]
  simp only [angleA, angleB, angleC]
  push_cast
  linear_combination hsum

omit [Fintype V] in
/-- **The index of every face is an integer** (in exact arithmetic), for a field without zeros. -/
theorem faceIndex_isInt (hM : M.Valid) (n : ℕ) {φ : F × Fin 3 → ℝ} (hφ : M.SignpostFrames φ)
    {u : V → ℂ} (hu : ∀ v, u v ≠ 0) (f : F) : ∃ m : ℤ, M.faceIndex n φ u f = m := by
  obtain ⟨m, hm⟩ := M.sum_transportAngle hM hφ f
  set r : Fin 3 → ℂ := fun k => Complex.exp (((n * M.transportAngle φ (f, k) : ℝ) : ℂ) * I)
  have hz : (u (M.b f) / (r 0 * u (M.a f))) * (u (M.c f) / (r 1 * u (M.b f))) *
      (u (M.a f) / (r 2 * u (M.c f))) *
        Complex.exp (((n * ∑ k, M.transportAngle φ (f, k) : ℝ) : ℂ) * I) = 1 := by
    have hr : Complex.exp (((n * ∑ k, M.transportAngle φ (f, k) : ℝ) : ℂ) * I) = r 0 * r 1 * r 2 := by
      simp only [r, ← Complex.exp_add, Fin.sum_univ_three]
      congr 1
      push_cast
      ring
    rw [hr]
    have := hu (M.a f)
    have := hu (M.b f)
    have := hu (M.c f)
    have hr0 : ∀ k, r k ≠ 0 := fun k => Complex.exp_ne_zero _
    have := hr0 0
    have := hr0 1
    have := hr0 2
    field_simp
  obtain ⟨m', hm'⟩ := arg_add_arg_add_arg_eq hz
  refine ⟨m' - n * m, ?_⟩
  unfold faceIndex
  rw [Fin.sum_univ_three]
  simp only [omega, r] at hm' ⊢
  simp only [halfedge_zero, halfedge_one, halfedge_two] at hm' ⊢
  rw [hm] at hm'
  have hπ : (2 * π) ≠ 0 := by positivity
  rw [div_eq_iff hπ, Int.cast_sub, Int.cast_mul, Int.cast_natCast]
  linear_combination hm'

omit [DecidableEq V] [Fintype V] [Fintype F] in
/-- The corner angle of `f` at `v`, as a sum over the corners of `f`. -/
theorem angleAt_eq_sum [DecidableEq V] (hM : M.Valid) (f : F) (v : V) :
    M.angleAt f v = ∑ k, if M.corner f k = v then M.cornerAngleAt f k else 0 := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  simp only [Fin.sum_univ_three, corner, cornerAngleAt, angleAt, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  by_cases h1 : v = M.a f <;> by_cases h2 : v = M.b f <;> by_cases h3 : v = M.c f <;>
    simp_all [eq_comm]

/-- **The holonomies of the faces add up to `2πχ`** on a closed valid manifold mesh, every
vertex in some face (`Θ_v ≠ 0`). -/
theorem sum_faceHolonomy (hM : M.Valid) (hman : M.Manifold) (hclosed : M.boundaryVertices = ∅)
    (hΘ : ∀ v, M.angleSum v ≠ 0) : ∑ f, M.faceHolonomy f = 2 * π * (M.eulerChar : ℝ) := by
  have h1 : ∑ f, M.faceHolonomy f = ∑ v, (M.tangentScale v - 1) * M.angleSum v := by
    simp only [angleSum, Finset.mul_sum, M.angleAt_eq_sum hM, mul_ite, mul_zero]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun f _ => ?_
    rw [Finset.sum_comm]
    simp only [faceHolonomy, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  rw [h1, ← M.gauss_bonnet_closed hM hman hclosed]
  refine Finset.sum_congr rfl fun v _ => ?_
  rw [tangentScale, angleDefect, Apps.scaled_minus_sum _ (hΘ v)]

/-- **The index sum, exactly** (`singularities` with `ω ∈ (−π, π]`): on a closed valid oriented
manifold mesh with signpost frames and a field without zeros,
`Σ_f index f = n χ + #{h : ω_h = π}/2`; each edge with `ω = π` (both of its halfedges, by
`omega_twin`) adds `+1` to the Poincaré–Hopf value `n χ`. -/
theorem sum_faceIndex (hM : M.Valid) (hman : M.Manifold) (hclosed : M.boundaryVertices = ∅)
    (hO : M.Oriented) (hΘ : ∀ v, M.angleSum v ≠ 0) (n : ℕ) {φ : F × Fin 3 → ℝ}
    {u : V → ℂ} (hu : ∀ v, u v ≠ 0) :
    ∑ f, M.faceIndex n φ u f =
      n * M.eulerChar + ((Finset.univ.filter fun h => M.omega n φ u h = π).card : ℝ) / 2 := by
  unfold faceIndex
  rw [← Finset.sum_div, Finset.sum_add_distrib, ← Finset.mul_sum, ← Fintype.sum_prod_type']
  set ω := M.omega n φ u
  have htw : ∀ h, ω (M.twinPerm hO h) = if ω h = π then π else -ω h := fun h =>
    M.omega_twin hO n φ hu (M.hasTwin_of_closed hM hO hclosed h)
  have hpair : 2 * ∑ h, ω h = 2 * π * (Finset.univ.filter fun h => ω h = π).card := by
    have hc : ∑ h, ω (M.twinPerm hO h) = ∑ h, ω h := Equiv.sum_comp _ ω
    have : ∑ h, (ω h + ω (M.twinPerm hO h)) = ∑ h, if ω h = π then 2 * π else 0 := by
      refine Finset.sum_congr rfl fun h _ => ?_
      rw [htw]
      split_ifs with e
      · rw [e]; ring
      · ring
    rw [Finset.sum_add_distrib, hc, Finset.sum_ite, Finset.sum_const_zero, add_zero,
      Finset.sum_const, nsmul_eq_mul] at this
    linarith
  have hΩ := M.sum_faceHolonomy hM hman hclosed hΘ
  rw [hΩ]
  simp only [Prod.mk.eta]
  have hπ : (2 * π) ≠ 0 := by positivity
  rw [div_eq_iff hπ]
  linear_combination hpair / 2

/-- **What `singularities` reported as `index_sum` before the fix**: the sum of the rounded indices. With
signpost frames every index is an integer (`faceIndex_isInt`), so rounding changes nothing and
`index_sum = n χ + #{h : ω_h = π}/2`. -/
theorem sum_round_faceIndex (hM : M.Valid) (hman : M.Manifold) (hclosed : M.boundaryVertices = ∅)
    (hO : M.Oriented) (hΘ : ∀ v, M.angleSum v ≠ 0) (n : ℕ) {φ : F × Fin 3 → ℝ}
    (hφ : M.SignpostFrames φ) {u : V → ℂ} (hu : ∀ v, u v ≠ 0) :
    ∑ f, ((round (M.faceIndex n φ u f) : ℤ) : ℝ) =
      n * M.eulerChar + ((Finset.univ.filter fun h => M.omega n φ u h = π).card : ℝ) / 2 := by
  rw [← M.sum_faceIndex hM hman hclosed hO hΘ n hu]
  refine Finset.sum_congr rfl fun f _ => ?_
  obtain ⟨m, hm⟩ := M.faceIndex_isInt hM n hφ hu f
  rw [hm, round_intCast]

/-- **Discrete Poincaré–Hopf** (GODF Eq. 19, App. B): if no halfedge has `ω = π` exactly, the
indices add up to `n χ`. -/
theorem sum_faceIndex_of_ne_pi (hM : M.Valid) (hman : M.Manifold)
    (hclosed : M.boundaryVertices = ∅) (hO : M.Oriented) (hΘ : ∀ v, M.angleSum v ≠ 0) (n : ℕ)
    {φ : F × Fin 3 → ℝ} {u : V → ℂ} (hu : ∀ v, u v ≠ 0) (hπ : ∀ h, M.omega n φ u h ≠ π) :
    ∑ f, M.faceIndex n φ u f = n * M.eulerChar := by
  rw [M.sum_faceIndex hM hman hclosed hO hΘ n hu]
  simp [hπ]

/-! ## The rule of the code: the twin gets `−ω`

`singularities` evaluates `ω` on the halfedge with the lower index and gives its twin `−ω`
(`Some(t) if t < h => -omega(t)`). `idx` is the numbering of the halfedges (`3f + k` in the
code; any injective numbering). -/

open Classical in
/-- `ω` as the code uses it: `−ω_{twin h}` if the twin has the lower index, else `ω_h`. -/
def omegaFixed (idx : F × Fin 3 → ℕ) (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ)
    (h : F × Fin 3) : ℝ :=
  if M.HasTwin h ∧ idx (M.twin h) < idx h then -M.omega n φ u (M.twin h) else M.omega n φ u h

/-- The index of `f` with `omegaFixed`. -/
def faceIndexFixed (idx : F × Fin 3 → ℕ) (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ) (f : F) : ℝ :=
  (∑ k, M.omegaFixed idx n φ u (f, k) + n * M.faceHolonomy f) / (2 * π)

omit [DecidableEq V] [Fintype V] [Fintype F] in
/-- The two halfedges of an edge get opposite angles (exactly, also at a tie). -/
theorem omegaFixed_add_twin (hM : M.Valid) (hO : M.Oriented) {idx : F × Fin 3 → ℕ}
    (hidx : Function.Injective idx) (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ) {h : F × Fin 3}
    (ht : M.HasTwin h) : M.omegaFixed idx n φ u h + M.omegaFixed idx n φ u (M.twin h) = 0 := by
  have htt := M.hasTwin_twin ht
  have hne : idx (M.twin h) ≠ idx h := fun e => M.twin_ne hM ht (hidx e)
  unfold omegaFixed
  rw [M.twin_twin hO]
  rcases lt_or_gt_of_ne hne with hl | hl
  · rw [ite_eq_left ⟨ht, hl⟩, ite_eq_right (fun hc => absurd hc.2 (lt_asymm hl))]
    ring
  · rw [ite_eq_right (fun hc => absurd hc.2 (lt_asymm hl)), ite_eq_left ⟨htt, hl⟩]
    ring

open Classical in
/-- `ω` literally as `singularities` evaluates it: `arg` on every halfedge; on the halfedge with
the higher index of a twin pair, shifted by `∓2π` when the pair sum `ω_h + ω_{twin h}` lies
outside `[−π, π]` (`if pair > PI { w - TAU } else if pair < -PI { w + TAU } else { w }`). -/
def omegaCode (idx : F × Fin 3 → ℕ) (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ)
    (h : F × Fin 3) : ℝ :=
  if M.HasTwin h ∧ idx (M.twin h) < idx h then
    (if π < M.omega n φ u h + M.omega n φ u (M.twin h) then M.omega n φ u h - 2 * π
      else if M.omega n φ u h + M.omega n φ u (M.twin h) < -π then M.omega n φ u h + 2 * π
      else M.omega n φ u h)
  else M.omega n φ u h

omit [DecidableEq V] [Fintype V] [Fintype F] in
/-- **The code's rule is `omegaFixed`** (in exact arithmetic): the shift happens exactly at a
tie, and there it turns `+π` into `−π = −ω_{twin h}`.

[formula: app-vector-singularity-index] -/
theorem omegaCode_eq_omegaFixed (hO : M.Oriented) (idx : F × Fin 3 → ℕ) (n : ℕ)
    (φ : F × Fin 3 → ℝ) {u : V → ℂ} (hu : ∀ v, u v ≠ 0) (h : F × Fin 3) :
    M.omegaCode idx n φ u h = M.omegaFixed idx n φ u h := by
  unfold omegaCode omegaFixed
  split_ifs with hc h1 h2
  · have e := M.omega_twin hO n φ hu hc.1
    split_ifs at e with hπ
    · rw [e, hπ]; ring
    · rw [e] at h1; linarith [Real.pi_pos]
  · have e := M.omega_twin hO n φ hu hc.1
    split_ifs at e with hπ
    · rw [e, hπ] at h2; linarith [Real.pi_pos]
    · rw [e] at h2; linarith [Real.pi_pos]
  · have e := M.omega_twin hO n φ hu hc.1
    split_ifs at e with hπ
    · rw [e, hπ] at h1; linarith [Real.pi_pos]
    · rw [e]; ring
  · rfl

omit [DecidableEq V] [Fintype V] [Fintype F] in
/-- `omegaFixed` differs from `ω` by full turns (by `−2π` exactly at a tie on the higher
halfedge). -/
theorem omegaFixed_eq (hO : M.Oriented) (idx : F × Fin 3 → ℕ) (n : ℕ) (φ : F × Fin 3 → ℝ)
    {u : V → ℂ} (hu : ∀ v, u v ≠ 0) (h : F × Fin 3) :
    ∃ m : ℤ, M.omegaFixed idx n φ u h = M.omega n φ u h + 2 * π * m := by
  unfold omegaFixed
  split_ifs with hc
  · have htt := M.hasTwin_twin hc.1
    have e := M.omega_twin hO n φ hu htt
    rw [M.twin_twin hO] at e
    -- `e : ω h = if ω (twin h) = π then π else −ω (twin h)`
    split_ifs at e with hπ
    · exact ⟨-1, by rw [hπ, e]; push_cast; ring⟩
    · exact ⟨0, by rw [e]; push_cast; ring⟩
  · exact ⟨0, by push_cast; ring⟩

omit [Fintype V] in
/-- With the code's rule every face index is still an integer. -/
theorem faceIndexFixed_isInt (hM : M.Valid) (hO : M.Oriented) (idx : F × Fin 3 → ℕ) (n : ℕ)
    {φ : F × Fin 3 → ℝ} (hφ : M.SignpostFrames φ) {u : V → ℂ} (hu : ∀ v, u v ≠ 0) (f : F) :
    ∃ m : ℤ, M.faceIndexFixed idx n φ u f = m := by
  obtain ⟨m, hm⟩ := M.faceIndex_isInt hM n hφ hu f
  obtain ⟨m0, h0⟩ := M.omegaFixed_eq hO idx n φ hu (f, 0)
  obtain ⟨m1, h1⟩ := M.omegaFixed_eq hO idx n φ hu (f, 1)
  obtain ⟨m2, h2⟩ := M.omegaFixed_eq hO idx n φ hu (f, 2)
  refine ⟨m + m0 + m1 + m2, ?_⟩
  unfold faceIndexFixed
  unfold faceIndex at hm
  rw [Fin.sum_univ_three] at hm ⊢
  rw [h0, h1, h2]
  have hπ : (2 * π) ≠ 0 := by positivity
  rw [div_eq_iff hπ] at hm ⊢
  push_cast
  linear_combination hm

omit [DecidableEq V] [Fintype V] [Fintype F] in
theorem cornerAngleAt_pos (hM : M.Valid) (f : F) (k : Fin 3) : 0 < M.cornerAngleAt f k := by
  fin_cases k
  · exact cornerAngle_pos (hM.tri f)
  · exact cornerAngle_pos (hM.tri f).rotate
  · exact cornerAngle_pos (hM.tri f).rotate.rotate

omit [Fintype V] in
/-- A vertex that is a corner of some face has a positive angle sum. -/
theorem angleSum_pos_of_corner (hM : M.Valid) {v : V} (hv : ∃ f k, M.corner f k = v) : 0 < M.angleSum v := by
  obtain ⟨f, k, hk⟩ := hv
  unfold angleSum
  simp only [M.angleAt_eq_sum hM]
  have hnn : ∀ g (l : Fin 3), 0 ≤ if M.corner g l = v then M.cornerAngleAt g l else 0 := by
    intro g l
    split_ifs
    · exact (M.cornerAngleAt_pos hM g l).le
    · exact le_rfl
  refine Finset.sum_pos' (fun g _ => Finset.sum_nonneg fun l _ => hnn g l) ⟨f, Finset.mem_univ _, ?_⟩
  refine Finset.sum_pos' (fun l _ => hnn f l) ⟨k, Finset.mem_univ _, ?_⟩
  rw [ite_eq_left hk]
  exact M.cornerAngleAt_pos hM f k

/-- **Poincaré–Hopf for the code's rule, unconditionally**: on a closed valid oriented manifold
mesh in which every vertex is a corner of some face, for a field without zeros and any
injective numbering of the halfedges, `Σ_f index f = n χ` — ties at `ω = π` included.

[formula: app-vector-singularity-index] -/
theorem sum_faceIndexFixed (hM : M.Valid) (hman : M.Manifold) (hclosed : M.boundaryVertices = ∅)
    (hO : M.Oriented) (hcorner : ∀ v, ∃ f k, M.corner f k = v) {idx : F × Fin 3 → ℕ}
    (hidx : Function.Injective idx) (n : ℕ) (φ : F × Fin 3 → ℝ) (u : V → ℂ) :
    ∑ f, M.faceIndexFixed idx n φ u f = n * M.eulerChar := by
  have hΘ : ∀ v, M.angleSum v ≠ 0 := fun v => (M.angleSum_pos_of_corner hM (hcorner v)).ne'
  unfold faceIndexFixed
  rw [← Finset.sum_div, Finset.sum_add_distrib, ← Finset.mul_sum, ← Fintype.sum_prod_type',
    M.sum_faceHolonomy hM hman hclosed hΘ]
  set ω := M.omegaFixed idx n φ u
  have hzero : ∑ h, ω h = 0 := by
    have hc : ∑ h, ω (M.twinPerm hO h) = ∑ h, ω h := Equiv.sum_comp _ ω
    have : ∑ h, (ω h + ω (M.twinPerm hO h)) = 0 :=
      Finset.sum_eq_zero fun h _ =>
        M.omegaFixed_add_twin hM hO hidx n φ u (M.hasTwin_of_closed hM hO hclosed h)
    rw [Finset.sum_add_distrib, hc] at this
    linarith
  simp only [Prod.mk.eta]
  rw [hzero, zero_add]
  have hπ : (2 * π) ≠ 0 := by positivity
  rw [div_eq_iff hπ]
  ring

/-- **What `singularities` reports as `index_sum`** (with the fix): the sum of the rounded
indices, which are integers for signpost frames (`faceIndexFixed_isInt`), is `n χ`.

[formula: app-vector-singularity-index] -/
theorem sum_round_faceIndexFixed (hM : M.Valid) (hman : M.Manifold)
    (hclosed : M.boundaryVertices = ∅) (hO : M.Oriented) (hcorner : ∀ v, ∃ f k, M.corner f k = v)
    {idx : F × Fin 3 → ℕ} (hidx : Function.Injective idx) (n : ℕ) {φ : F × Fin 3 → ℝ}
    (hφ : M.SignpostFrames φ) {u : V → ℂ} (hu : ∀ v, u v ≠ 0) :
    ∑ f, ((round (M.faceIndexFixed idx n φ u f) : ℤ) : ℝ) = n * M.eulerChar := by
  rw [← M.sum_faceIndexFixed hM hman hclosed hO hcorner hidx n φ u]
  refine Finset.sum_congr rfl fun f _ => ?_
  obtain ⟨m, hm⟩ := M.faceIndexFixed_isInt hM hO idx n hφ hu f
  rw [hm, round_intCast]

end TriMesh

end DDG
