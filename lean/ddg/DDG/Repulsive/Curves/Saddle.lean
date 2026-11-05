import DDG.Repulsive.Curves.Gram

/-!
# Repulsive Curves: the saddle systems, the projector and the scaling of the direction
(paper Sec. 5.3, Eq. 22, 23, 29, 30, 33; SPEC §6, D4, D5)

Generic in the index types: `A : Matrix n n ℝ` is the metric (`Ā = A ⊗ I₃` on the interleaved
coordinates, `Gram.lean`), `C : Matrix m n ℝ` the constraint Jacobian. `SaddleSol A C b r x μ`
means `[A Cᵀ; C 0][x; μ] = [b; r]`.

* `rc-saddle-gradient` (Eq. 22): for symmetric `A`, a solution `g̃` of the system with right side
  `(dFᵀ, 0)` satisfies `Φ(y) = Φ(g̃) + ½ (y − g̃)ᵀ A (y − g̃)` on `ker C`, `Φ(y) = ½ yᵀAy − dF·y`;
  equivalently (if `A g = dFᵀ`) `½‖y − g‖²_A = ½‖g̃ − g‖²_A + ½‖y − g̃‖²_A`: `g̃` is the `A`-closest
  point to `g` in `ker C`. If `A` is positive definite on `ker C` and `C` has full row rank, the
  solution exists and is unique (`saddle_exists`, `saddle_unique`); positive definiteness on `ker C` follows
  from `A` PSD with `ker A ∩ ker C = 0` (`posDefOn_of_ker`), which is the translation check of
  §5.2 when `ker Ā` is the translations (D5).
* `rc-descent-direction`: `dF·g̃ = g̃ᵀ A g̃`, positive unless `g̃ = 0`.
* `rc-chord-projection` (Eq. 23): a solution `x` of the system with right side `(0, −Φ)` has
  `yᵀAy = xᵀAx + (y − x)ᵀA(y − x)` for every `y` with `C y = −Φ`: the `A`-smallest solution of the
  linearized constraints.
* `rc-projector` (D4): `C† = Cᵀ(CCᵀ)⁻¹`, `P = I − C†C` for `CCᵀ` invertible: `C C† = I` (so the
  paper's `P := CC†` is the identity), `C P = 0`, `P² = P`, `Pᵀ = P`, `range P = ker C`.
* `rc-projected-system` (Eq. 30): `PᵀAPy = Pᵀa` iff `x = P y` solves Eq. 29 for some `λ`.
* `rc-projection-transform` (Eq. 33): `z = C†b`, `(y, µ)` solving Eq. 29 with right side `A z`
  ⟹ `(z − y, −µ)` solves Eq. 23 with lower block `b`.
* `rc-mean-free-direction` (§6.4): the rows `M_c = [m₁I₃ … m_|V|I₃]`: `M_c z = 0` iff `z` is
  `L²`-orthogonal to the constant fields; `M_c M_cᵀ = (Σ m_v²) I₃`; `M_c t = (Σ m_v) t` on a
  translation, so for `Σ m_v ≠ 0` and `ker A` = constants the system is uniquely solvable.
* `rc-gradient-scaling` (§4.5): if `A ↦ sA`, `dF ↦ t dF` and `C` is unchanged, the direction
  becomes `(t/s) g̃`; with `s = c^{−(2σ+1)}` (`gramB_smul`, `gramB0_smul`) and `t = c^{1+α−β}`
  (the differential of a function homogeneous of degree `2 + α − β`, `grad_homogeneous`), the
  direction scales by `c^{2σ+2+α−β}`.
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Finset
open scoped Kronecker

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

/-- `(x, μ)` solves `[A Cᵀ; C 0][x; μ] = [b; r]`. -/
def SaddleSol (A : Matrix n n ℝ) (C : Matrix m n ℝ) (b : n → ℝ) (r : m → ℝ) (x : n → ℝ)
    (μ : m → ℝ) : Prop :=
  A *ᵥ x + Cᵀ *ᵥ μ = b ∧ C *ᵥ x = r

/-- `A` is positive definite on `ker C`. -/
def PosDefOn (A : Matrix n n ℝ) (C : Matrix m n ℝ) : Prop :=
  ∀ x, C *ᵥ x = 0 → x ≠ 0 → 0 < x ⬝ᵥ A *ᵥ x

section basic

variable {A : Matrix n n ℝ} {C : Matrix m n ℝ}

omit [DecidableEq n] [DecidableEq m] in
theorem transpose_dot (C : Matrix m n ℝ) (μ : m → ℝ) (y : n → ℝ) :
    (Cᵀ *ᵥ μ) ⬝ᵥ y = μ ⬝ᵥ C *ᵥ y := by
  rw [mulVec_transpose, dotProduct_mulVec]

omit [Fintype m] [DecidableEq n] [DecidableEq m] in
theorem sym_dot (hA : Aᵀ = A) (x y : n → ℝ) : x ⬝ᵥ A *ᵥ y = y ⬝ᵥ A *ᵥ x := by
  rw [dotProduct_mulVec, ← mulVec_transpose, hA, dotProduct_comm]

omit [Fintype m] [DecidableEq n] [DecidableEq m] in
theorem mulVec_dot (hA : Aᵀ = A) (x y : n → ℝ) : (A *ᵥ x) ⬝ᵥ y = x ⬝ᵥ A *ᵥ y := by
  rw [dotProduct_comm, sym_dot hA]

omit [DecidableEq n] [DecidableEq m] in
/-- **The saddle system of Eq. 22 characterizes the `A`-closest tangent direction**: for symmetric
`A` and a solution `(g̃, λ)` with right side `(b, 0)`, every `y ∈ ker C` has
`Φ(y) = Φ(g̃) + ½ (y − g̃)ᵀ A (y − g̃)` with `Φ(y) = ½ yᵀAy − b·y`.

[formula: rc-saddle-gradient] -/
theorem saddle_split (hA : Aᵀ = A) {b : n → ℝ} {g : n → ℝ} {μ : m → ℝ}
    (h : SaddleSol A C b 0 g μ) {y : n → ℝ} (hy : C *ᵥ y = 0) :
    (1 / 2) * (y ⬝ᵥ A *ᵥ y) - b ⬝ᵥ y =
      (1 / 2) * (g ⬝ᵥ A *ᵥ g) - b ⬝ᵥ g + (1 / 2) * ((y - g) ⬝ᵥ A *ᵥ (y - g)) := by
  obtain ⟨h1, h2⟩ := h
  have hby : b ⬝ᵥ y = g ⬝ᵥ A *ᵥ y := by
    rw [← h1, add_dotProduct, transpose_dot, hy, dotProduct_zero, add_zero, mulVec_dot hA]
  have hbg : b ⬝ᵥ g = g ⬝ᵥ A *ᵥ g := by
    rw [← h1, add_dotProduct, transpose_dot, h2, dotProduct_zero, add_zero, mulVec_dot hA]
  have hd : (y - g) ⬝ᵥ A *ᵥ (y - g) =
      y ⬝ᵥ A *ᵥ y - 2 * (g ⬝ᵥ A *ᵥ y) + g ⬝ᵥ A *ᵥ g := by
    rw [mulVec_sub, sub_dotProduct, dotProduct_sub, dotProduct_sub, sym_dot hA y g]
    ring
  rw [hby, hbg, hd]
  ring

omit [DecidableEq n] [DecidableEq m] in
/-- **`g̃` is the `A`-closest point of `ker C` to the unconstrained gradient `g`** (`A g = b`):
`½‖y − g‖²_A = ½‖g̃ − g‖²_A + ½‖y − g̃‖²_A` for every `y ∈ ker C`.

[formula: rc-saddle-gradient] -/
theorem saddle_closest (hA : Aᵀ = A) {b g0 g : n → ℝ} {μ : m → ℝ} (hg0 : A *ᵥ g0 = b)
    (h : SaddleSol A C b 0 g μ) {y : n → ℝ} (hy : C *ᵥ y = 0) :
    (1 / 2) * ((y - g0) ⬝ᵥ A *ᵥ (y - g0)) =
      (1 / 2) * ((g - g0) ⬝ᵥ A *ᵥ (g - g0)) + (1 / 2) * ((y - g) ⬝ᵥ A *ᵥ (y - g)) := by
  have hs := saddle_split hA h hy
  have e : ∀ z, (z - g0) ⬝ᵥ A *ᵥ (z - g0) = z ⬝ᵥ A *ᵥ z - 2 * (b ⬝ᵥ z) + g0 ⬝ᵥ A *ᵥ g0 := by
    intro z
    rw [mulVec_sub, sub_dotProduct, dotProduct_sub, dotProduct_sub, ← hg0, sym_dot hA z g0,
      mulVec_dot hA]
    ring
  rw [e y, e g]
  linarith

omit [DecidableEq n] [DecidableEq m] in
/-- For `A` positive definite on `ker C`, the solution `g̃` is the **unique minimizer** of `Φ` on
`ker C`.

[formula: rc-saddle-gradient] -/
theorem saddle_isMin (hA : Aᵀ = A) (hpd : PosDefOn A C) {b g : n → ℝ} {μ : m → ℝ}
    (h : SaddleSol A C b 0 g μ) {y : n → ℝ} (hy : C *ᵥ y = 0) (hne : y ≠ g) :
    (1 / 2) * (g ⬝ᵥ A *ᵥ g) - b ⬝ᵥ g < (1 / 2) * (y ⬝ᵥ A *ᵥ y) - b ⬝ᵥ y := by
  rw [saddle_split hA h hy]
  have hd : C *ᵥ (y - g) = 0 := by rw [mulVec_sub, hy, h.2, sub_zero]
  have := hpd (y - g) hd (sub_ne_zero.mpr hne)
  linarith

omit [DecidableEq n] [DecidableEq m] in
/-- **Uniqueness**: for `A` positive definite on `ker C`, two solutions of the same saddle system
have the same `x`; if moreover `Cᵀ` is injective (full row rank), the same multipliers.

[formula: rc-saddle-gradient] -/
theorem saddle_unique (hpd : PosDefOn A C) {b : n → ℝ} {r : m → ℝ} {x x' : n → ℝ}
    {μ μ' : m → ℝ} (h : SaddleSol A C b r x μ) (h' : SaddleSol A C b r x' μ') :
    x = x' ∧ ((∀ ν, Cᵀ *ᵥ ν = 0 → ν = 0) → μ = μ') := by
  have hd1 : A *ᵥ (x - x') + Cᵀ *ᵥ (μ - μ') = 0 := by
    rw [mulVec_sub, mulVec_sub, sub_add_sub_comm, h.1, h'.1, sub_self]
  have hd2 : C *ᵥ (x - x') = 0 := by rw [mulVec_sub, h.2, h'.2, sub_self]
  have hx : x - x' = 0 := by
    by_contra hne
    have hpos := hpd _ hd2 hne
    have : (x - x') ⬝ᵥ A *ᵥ (x - x') = 0 := by
      have := congrArg (fun z => (x - x') ⬝ᵥ z) hd1
      simp only [dotProduct_add, dotProduct_zero] at this
      rw [dotProduct_comm (x - x') (Cᵀ *ᵥ _), transpose_dot, hd2, dotProduct_zero,
        add_zero] at this
      exact this
    linarith
  refine ⟨sub_eq_zero.mp hx, fun hinj => ?_⟩
  rw [hx, mulVec_zero, zero_add] at hd1
  exact sub_eq_zero.mp (hinj _ hd1)

omit [DecidableEq m] in
/-- **Existence**: for `A` positive definite on `ker C` and `C` of full row rank (`Cᵀ`
injective), the saddle system has a solution for every right side (the block matrix is
invertible).

[formula: rc-saddle-gradient] -/
theorem saddle_exists [DecidableEq m] (hpd : PosDefOn A C) (hinj : ∀ ν, Cᵀ *ᵥ ν = 0 → ν = 0)
    (b : n → ℝ) (r : m → ℝ) : ∃ x μ, SaddleSol A C b r x μ := by
  set M : Matrix (n ⊕ m) (n ⊕ m) ℝ := fromBlocks A Cᵀ C 0 with hM
  have hsol : ∀ (b : n → ℝ) (r : m → ℝ) (v : n ⊕ m → ℝ), M *ᵥ v = Sum.elim b r ↔
      SaddleSol A C b r (v ∘ Sum.inl) (v ∘ Sum.inr) := by
    intro b r v
    rw [hM, fromBlocks_mulVec]
    constructor
    · intro h
      refine ⟨funext fun i => ?_, funext fun j => ?_⟩
      · have := congrFun h (Sum.inl i)
        simpa using this
      · have := congrFun h (Sum.inr j)
        simpa using this
    · rintro ⟨h1, h2⟩
      funext k
      cases k with
      | inl i => simpa using congrFun h1 i
      | inr j => simpa using congrFun h2 j
  have hinjM : Function.Injective M.mulVec := by
    intro v w hvw
    have hv : M *ᵥ (v - w) = 0 := by rw [mulVec_sub]; exact sub_eq_zero.mpr hvw
    have h0 : SaddleSol A C 0 0 ((v - w) ∘ Sum.inl) ((v - w) ∘ Sum.inr) := by
      rw [← hsol 0 0, hv]
      funext k; cases k <;> rfl
    have hz : SaddleSol A C 0 0 (0 : n → ℝ) (0 : m → ℝ) := by
      simp [SaddleSol]
    obtain ⟨hx, hμ⟩ := saddle_unique hpd h0 hz
    have hμ' := hμ hinj
    have : v - w = 0 := by
      funext k
      cases k with
      | inl i => exact congrFun hx i
      | inr j => exact congrFun hμ' j
    exact sub_eq_zero.mp this
  have hsurj : Function.Surjective M.mulVec :=
    mulVec_surjective_iff_isUnit.mpr (mulVec_injective_iff_isUnit.mp hinjM)
  obtain ⟨v, hv⟩ := hsurj (Sum.elim b r)
  exact ⟨_, _, (hsol b r v).mp hv⟩

omit [Fintype m] [DecidableEq n] [DecidableEq m] in
/-- **Positive definiteness on `ker C`** from a PSD `A` with `ker A ∩ ker C = 0` (e.g. `ker Ā` the
translations and `C` of full rank on them: the check `TranslationNotFixed` of SPEC §5.2, D5).

[formula: rc-saddle-gradient] -/
theorem posDefOn_of_ker (hA : A.PosSemidef) (hker : ∀ x, A *ᵥ x = 0 → C *ᵥ x = 0 → x = 0) :
    PosDefOn A C := by
  intro x hx hne
  rcases (psd_form_nonneg hA x).eq_or_lt with h | h
  · exact absurd (hker x ((mulVec_eq_zero_iff_form hA x).mpr h.symm) hx) hne
  · exact h

omit [DecidableEq n] [DecidableEq m] in
/-- **The projected gradient is a descent direction**: `dF·g̃ = g̃ᵀ A g̃`, and for `A` positive
definite on `ker C`, `dF·g̃ > 0` unless `g̃ = 0` (so `−g̃` decreases `F` and is tangent to the
constraint set).

[formula: rc-descent-direction] -/
theorem descent_direction (hA : Aᵀ = A) {b g : n → ℝ} {μ : m → ℝ} (h : SaddleSol A C b 0 g μ) :
    b ⬝ᵥ g = g ⬝ᵥ A *ᵥ g ∧ (PosDefOn A C → g ≠ 0 → 0 < b ⬝ᵥ g) := by
  have e : b ⬝ᵥ g = g ⬝ᵥ A *ᵥ g := by
    rw [← h.1, add_dotProduct, transpose_dot, h.2, dotProduct_zero, add_zero, mulVec_dot hA]
  exact ⟨e, fun hpd hne => e ▸ hpd g h.2 hne⟩

omit [DecidableEq n] [DecidableEq m] in
/-- **One chord iteration is the `A`-least-norm solution of the linearized constraints** (Eq. 23):
if `(x, μ)` solves the system with right side `(0, r)` (`r = −Φ(γ̃)`), then every `y` with
`C y = r` has `yᵀAy = xᵀAx + (y − x)ᵀA(y − x)`, so `xᵀAx ≤ yᵀAy` for `A` PSD on `ker C`.

[formula: rc-chord-projection] -/
theorem chord_least_norm (hA : Aᵀ = A) {r : m → ℝ} {x : n → ℝ} {μ : m → ℝ}
    (h : SaddleSol A C 0 r x μ) {y : n → ℝ} (hy : C *ᵥ y = r) :
    y ⬝ᵥ A *ᵥ y = x ⬝ᵥ A *ᵥ x + (y - x) ⬝ᵥ A *ᵥ (y - x) := by
  have hd : C *ᵥ (y - x) = 0 := by rw [mulVec_sub, hy, h.2, sub_self]
  have hAx : A *ᵥ x = -(Cᵀ *ᵥ μ) := eq_neg_of_add_eq_zero_left h.1
  have hx : x ⬝ᵥ A *ᵥ (y - x) = 0 := by
    rw [← mulVec_dot hA, hAx, neg_dotProduct, transpose_dot, hd, dotProduct_zero, neg_zero]
  have e : y = x + (y - x) := by abel
  conv_lhs => rw [e]
  rw [mulVec_add, add_dotProduct, dotProduct_add, dotProduct_add, hx, sym_dot hA (y - x) x, hx]
  ring

end basic

/-! ## The null-space projector (D4) -/

section projector

variable (C : Matrix m n ℝ)

/-- `C† = Cᵀ (C Cᵀ)⁻¹`. -/
def pinv : Matrix n m ℝ := Cᵀ * (C * Cᵀ)⁻¹

/-- `P = I − C† C`.

[formula: rc-projector] -/
def projP : Matrix n n ℝ := 1 - pinv C * C

variable {C}

/-- **The projector of D4** (for `C Cᵀ` invertible, i.e. `C` of full row rank): `C C† = I` (so
the printed `P := C C†` is the identity), `C P = 0`, `P² = P`, `Pᵀ = P`, and `range P = ker C`.

[formula: rc-projector] -/
theorem projector_props (hC : IsUnit (C * Cᵀ).det) :
    C * pinv C = 1 ∧ C * projP C = 0 ∧ projP C * projP C = projP C ∧ (projP C)ᵀ = projP C ∧
      (∀ x, C *ᵥ (projP C *ᵥ x) = 0) ∧ (∀ x, C *ᵥ x = 0 → projP C *ᵥ x = x) := by
  have h1 : C * pinv C = 1 := by
    rw [pinv, ← Matrix.mul_assoc, mul_nonsing_inv _ hC]
  have h2 : C * projP C = 0 := by
    rw [projP, Matrix.mul_sub, Matrix.mul_one, ← Matrix.mul_assoc, h1, Matrix.one_mul, sub_self]
  have hQ : pinv C * C * (pinv C * C) = pinv C * C := by
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc C, h1, Matrix.one_mul]
  refine ⟨h1, h2, ?_, ?_, fun x => ?_, fun x hx => ?_⟩
  · rw [projP, Matrix.sub_mul, Matrix.mul_sub, Matrix.mul_sub, Matrix.one_mul, Matrix.mul_one,
      Matrix.one_mul, hQ]
    abel
  · have hs : (C * Cᵀ)ᵀ = C * Cᵀ := by rw [transpose_mul, transpose_transpose]
    rw [projP, transpose_sub, transpose_one, pinv, transpose_mul, transpose_mul,
      transpose_nonsing_inv, hs, transpose_transpose, Matrix.mul_assoc]
  · rw [mulVec_mulVec, h2, zero_mulVec]
  · rw [projP, sub_mulVec, one_mulVec, ← mulVec_mulVec, hx, mulVec_zero, sub_zero]

/-- **Eq. 30 ⟺ Eq. 29**: `PᵀAPy = Pᵀa` iff `x = P y` solves `[A Cᵀ; C 0][x; λ] = [a; 0]` for some
`λ` (and `x` is unique for `A` positive definite on `ker C`, `saddle_unique`).

[formula: rc-projected-system] -/
theorem projected_system (hC : IsUnit (C * Cᵀ).det) (A : Matrix n n ℝ) (a y : n → ℝ) :
    (projP C)ᵀ *ᵥ (A *ᵥ (projP C *ᵥ y)) = (projP C)ᵀ *ᵥ a ↔
      ∃ μ, SaddleSol A C a 0 (projP C *ᵥ y) μ := by
  obtain ⟨h1, h2, -, hsym, hran, hfix⟩ := projector_props hC
  have hPC : (projP C)ᵀ * Cᵀ = 0 := by rw [← transpose_mul, h2, transpose_zero]
  constructor
  · intro h
    set r := a - A *ᵥ (projP C *ᵥ y) with hr
    have hPr' : (projP C)ᵀ *ᵥ r = 0 := by rw [hr, mulVec_sub, h, sub_self]
    have hPr : projP C *ᵥ r = 0 := by rwa [hsym] at hPr'
    -- `r = C† C r`
    have hr' : r = pinv C *ᵥ (C *ᵥ r) := by
      have : projP C *ᵥ r = r - pinv C *ᵥ (C *ᵥ r) := by
        rw [projP, sub_mulVec, one_mulVec, mulVec_mulVec]
      rw [hPr] at this
      exact (sub_eq_zero.mp this.symm)
    refine ⟨(C * Cᵀ)⁻¹ *ᵥ (C *ᵥ r), ?_, hran y⟩
    have e : Cᵀ *ᵥ ((C * Cᵀ)⁻¹ *ᵥ (C *ᵥ r)) = r := by
      rw [mulVec_mulVec, ← pinv]; exact hr'.symm
    rw [e, hr]
    abel
  · rintro ⟨μ, hs1, -⟩
    have e : (projP C)ᵀ *ᵥ (Cᵀ *ᵥ μ) = 0 := by rw [mulVec_mulVec, hPC, zero_mulVec]
    rw [← hs1, mulVec_add, e, add_zero]

/-- **Eq. 33**: with `z = C† b` and `(y, µ)` solving `[A Cᵀ; C 0][y; µ] = [A z; 0]`, the vector
`x = z − y` (with multipliers `−µ`) solves Eq. 23, `[A Cᵀ; C 0][x; ·] = [0; b]`.

[formula: rc-projection-transform] -/
theorem projection_transform (hC : IsUnit (C * Cᵀ).det) (A : Matrix n n ℝ) (b : m → ℝ)
    {y : n → ℝ} {μ : m → ℝ} (h : SaddleSol A C (A *ᵥ (pinv C *ᵥ b)) 0 y μ) :
    SaddleSol A C 0 b (pinv C *ᵥ b - y) (-μ) := by
  obtain ⟨h1, -⟩ := projector_props hC
  constructor
  · rw [mulVec_sub, mulVec_neg, ← h.1]
    abel
  · rw [mulVec_sub, h.2, sub_zero, mulVec_mulVec, h1, one_mulVec]

end projector

/-! ## The rows of the soft-constraint methods (§6.4) -/

section meanFree

variable {V : Type*} [Fintype V] [DecidableEq V] (mass : V → ℝ)

/-- `M_c = [m₁ I₃ … m_|V| I₃]` on the interleaved coordinates `(v, c)`. -/
def massRows : Matrix (Fin 3) (V × Fin 3) ℝ := fun c vc => if c = vc.2 then mass vc.1 else 0

omit [DecidableEq V] in
theorem massRows_mulVec (z : V × Fin 3 → ℝ) (c : Fin 3) :
    (massRows mass *ᵥ z) c = ∑ v, mass v * z (v, c) := by
  simp only [massRows, mulVec, dotProduct, Fintype.sum_prod_type, ite_mul, zero_mul,
    Finset.sum_ite_eq, Finset.mem_univ, ite_true]

omit [DecidableEq V] in
/-- **The rows `M_c` of §6.4**: `M_c z = 0` iff `z` is `L²`-orthogonal (`⟨z, t⟩ = Σ_v m_v z_v·t`) to
every constant field `t`; `M_c M_cᵀ = (Σ_v m_v²) I₃` (full rank for nonzero masses); a translation
`t` has `M_c t = (Σ_v m_v) t`.

[formula: rc-mean-free-direction] -/
theorem massRows_props (z : V × Fin 3 → ℝ) (t : V3) :
    (massRows mass *ᵥ z = 0 ↔ ∀ s : V3, ∑ v, mass v * ∑ c, z (v, c) * s c = 0) ∧
      massRows mass * (massRows mass)ᵀ = (∑ v, mass v ^ 2) • (1 : Matrix (Fin 3) (Fin 3) ℝ) ∧
      massRows mass *ᵥ (fun vc => t vc.2) = (∑ v, mass v) • t := by
  refine ⟨⟨fun h s => ?_, fun h => ?_⟩, ?_, ?_⟩
  · have hc : ∀ c, ∑ v, mass v * z (v, c) = 0 := fun c => by
      rw [← massRows_mulVec, h]; rfl
    simp only [Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_eq_zero fun c _ => ?_
    calc ∑ v, mass v * (z (v, c) * s c) = (∑ v, mass v * z (v, c)) * s c := by
          rw [Finset.sum_mul]; exact Finset.sum_congr rfl fun v _ => by ring
      _ = 0 := by rw [hc c, zero_mul]
  · funext c
    rw [massRows_mulVec]
    have := h (Pi.single c 1)
    simpa [Pi.single_apply, mul_comm] using this
  · ext c c'
    simp only [mul_apply, transpose_apply, massRows, Fintype.sum_prod_type, Matrix.smul_apply, one_apply,
      smul_eq_mul, mul_ite, ite_mul, zero_mul, mul_zero]
    by_cases h : c = c'
    · subst h
      simp [sq]
    · simp [h]
  · funext c
    rw [massRows_mulVec, Pi.smul_apply, smul_eq_mul, Finset.sum_mul]

/-- **The system of §6.4 is uniquely solvable**: if `ker A` is the constants and `Σ_v m_v ≠ 0`,
then `Ā = A ⊗ I₃` (PSD) is positive definite on `ker M_c`.

[formula: rc-mean-free-direction] -/
theorem massRows_posDefOn {A : Matrix V V ℝ} (hA : A.PosSemidef)
    (hker : ∀ u : V → ℝ, A *ᵥ u = 0 ↔ ∃ c, ∀ v, u v = c) (hm : ∑ v, mass v ≠ 0) :
    PosDefOn (A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) (massRows mass) := by
  have hpsd : (A ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)).PosSemidef := by
    refine posSemidef_of_form (fun u v => ?_) (fun u => ?_)
    · rw [dot_kron_one, dot_kron_one, Finset.sum_comm]
      refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
      have hs : A J I = A I J := by
        have := hA.1
        rw [IsHermitian, conjTranspose_eq_transpose_of_trivial] at this
        rw [← transpose_apply A I J, this]
      rw [hs]
      congr 1
      exact Finset.sum_congr rfl fun c _ => by ring
    · rw [dot_kron_one]
      have : ∑ I, ∑ J, A I J * ∑ c, u (I, c) * u (J, c) =
          ∑ c, (fun v => u (v, c)) ⬝ᵥ A *ᵥ (fun v => u (v, c)) := by
        have e1 : ∀ c, (fun v => u (v, c)) ⬝ᵥ A *ᵥ (fun v => u (v, c)) =
            ∑ I, ∑ J, A I J * (u (I, c) * u (J, c)) := by
          intro c
          simp only [dotProduct, mulVec, Finset.mul_sum]
          exact Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => by ring
        simp only [e1, Finset.mul_sum]
        symm
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun I _ => ?_
        rw [Finset.sum_comm]
      rw [this]
      exact Finset.sum_nonneg fun c _ => psd_form_nonneg hA _
  refine posDefOn_of_ker hpsd fun x hx hc => ?_
  obtain ⟨t, ht⟩ := (kron_ker_translations A hker x).mp hx
  have hx' : x = fun vc => t vc.2 := funext fun vc => ht vc.1 vc.2
  rw [hx', (massRows_props mass x t).2.2] at hc
  have ht0 : t = 0 := by
    funext c
    have := congrFun hc c
    simp only [Pi.smul_apply, smul_eq_mul, Pi.zero_apply, mul_eq_zero] at this
    exact this.resolve_left hm
  rw [hx', ht0]
  rfl

end meanFree

/-! ## The scaling of the direction (§4.5) -/

section scaling

variable {A : Matrix n n ℝ} {C : Matrix m n ℝ}

omit [DecidableEq n] [DecidableEq m] in
/-- Scaling the metric by `s ≠ 0` and the differential by `t` scales the direction by `t/s`
(the constraint rows unchanged). -/
theorem saddle_scale {b g : n → ℝ} {μ : m → ℝ} (h : SaddleSol A C b 0 g μ) {s : ℝ} (hs : s ≠ 0)
    (t : ℝ) : SaddleSol (s • A) C (t • b) 0 ((t / s) • g) (t • μ) := by
  constructor
  · rw [smul_mulVec, mulVec_smul, mulVec_smul, smul_smul, ← h.1, smul_add]
    congr 1
    rw [mul_div_cancel₀ _ hs]
  · rw [mulVec_smul, h.2, smul_zero]

omit [DecidableEq n] [DecidableEq m] in
/-- **The direction scales by `c^{2σ+2+α−β}`** (SPEC §4.5; A.4.2 corrected by D3): if the metric
of `cγ` is `c^{−(2σ+1)}` times that of `γ` (`gramB_smul`, `gramB0_smul`), the differential is
`c^{1+α−β}` times (`grad_homogeneous` with degree `2 + α − β`), the Jacobian rows are unchanged
(rows homogeneous of degree 0), and `A` is positive definite on `ker C`, then the projected
direction of `cγ` is `c^{2σ+2+α−β}` times that of `γ`.

[formula: rc-gradient-scaling] -/
theorem direction_scaling (hpd : PosDefOn A C) {σ α β c : ℝ} (hc : 0 < c) {b g g' : n → ℝ}
    {μ μ' : m → ℝ} (h : SaddleSol A C b 0 g μ)
    (h' : SaddleSol (c ^ (-(2 * σ + 1)) • A) C (c ^ (1 + α - β) • b) 0 g' μ') :
    g' = c ^ (2 * σ + 2 + α - β) • g := by
  have hs : c ^ (-(2 * σ + 1)) ≠ 0 := (Real.rpow_pos_of_pos hc _).ne'
  have hpd' : PosDefOn (c ^ (-(2 * σ + 1)) • A) C := by
    intro x hx hne
    rw [smul_mulVec, dotProduct_smul, smul_eq_mul]
    exact mul_pos (Real.rpow_pos_of_pos hc _) (hpd x hx hne)
  have := (saddle_unique hpd' h' (saddle_scale h hs (c ^ (1 + α - β)))).1
  rw [this, ← Real.rpow_sub hc]
  congr 2
  ring

/-- **The differential of a homogeneous function**: if `F(cγ) = c^p F(γ)` for every `γ` and the
derivatives of `F` along curves through `γ` and through `cγ` are given by `g` and `g'`, then
`g' = c^{p−1} g` (for `F = Ê`, `p = 2 + α − β`, `rc-energy-scaling`: `dÊ(cγ) = c^{1+α−β} dÊ(γ)`).

[formula: rc-gradient-scaling] -/
theorem grad_homogeneous {V : Type*} [Fintype V] [DecidableEq V] {F : (V → V3) → ℝ} {p c : ℝ} (hc : 0 < c)
    (hhom : ∀ γ, F (c • γ) = c ^ p * F γ) {γ g g' : V → V3}
    (hd : ∀ (Γ : ℝ → V → V3) (Γ' : V → V3), Γ 0 = γ →
      (∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) 0) →
        HasDerivAt (fun s => F (Γ s)) (∑ v, g v ⬝ᵥ Γ' v) 0)
    (hd' : ∀ (Γ : ℝ → V → V3) (Γ' : V → V3), Γ 0 = c • γ →
      (∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) 0) →
        HasDerivAt (fun s => F (Γ s)) (∑ v, g' v ⬝ᵥ Γ' v) 0) :
    g' = c ^ (p - 1) • g := by
  have key : ∀ δ : V → V3, ∑ v, g' v ⬝ᵥ δ v = c ^ (p - 1) * ∑ v, g v ⬝ᵥ δ v := by
    intro δ
    have hcurve : ∀ v, HasDerivAt (fun s : ℝ => γ v + s • δ v) (δ v) 0 := fun v => by
      simpa using ((hasDerivAt_id (0 : ℝ)).smul_const (δ v)).const_add (γ v)
    have h1 := hd (fun s v => γ v + s • δ v) δ (by funext v; simp) hcurve
    have h2 := hd' (fun s v => c • (γ v + s • δ v)) (fun v => c • δ v) (by funext v; simp)
      (fun v => (hcurve v).const_smul c)
    have h3 : HasDerivAt (fun s : ℝ => F (fun v => c • (γ v + s • δ v)))
        (c ^ p * ∑ v, g v ⬝ᵥ δ v) 0 := by
      have := h1.const_mul (c ^ p)
      refine this.congr_of_eventuallyEq (Filter.Eventually.of_forall fun s => ?_)
      have : (fun v => c • (γ v + s • δ v)) = c • fun v => γ v + s • δ v := rfl
      simp only [this, hhom]
    have := h2.unique h3
    simp only [dotProduct_smul, smul_eq_mul, ← Finset.mul_sum] at this
    have hcp : c ^ p = c * c ^ (p - 1) := by
      rw [show p = 1 + (p - 1) by ring, Real.rpow_add hc, Real.rpow_one]
      ring_nf
    rw [hcp, mul_assoc] at this
    exact mul_left_cancel₀ hc.ne' this
  funext v
  funext i
  have := key (fun w => if w = v then Pi.single i 1 else 0)
  simp only [apply_ite, dotProduct_zero, Finset.sum_ite_eq', Finset.mem_univ,
    ite_true] at this
  simpa [dotProduct, Pi.single_apply] using this

end scaling

end RC

end DDG
