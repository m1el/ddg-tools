import DdgProofs.Developability.Vec

/-!
# Developability: the covariance energy `E_λ` (paper Eq. 3, 4, 10)

`crates/ddg-developability/SPEC.md` §3.4, App. B.3 of the paper. The oracle's `Ddg.Dev.Sym3`,
`covMatrix` (`A_v = Σ θ_σ N_σ N_σᵀ`), `covObjective` (`Σ θ_σ ⟨u, N_σ⟩²`), `covSensitivityN`,
`covSensitivityTheta` at `ℝ`.

* `quad_covMatrix` — `uᵀ A_v u = Σ θ ⟨u, N⟩²` (Eq. 4 is the matrix of the quadratic form of Eq. 3);
* `covObjective_isLeast` — if `e₀, e₁, e₂` is an orthonormal eigenbasis of `A_v` with eigenvalues
  `λ₀ ≤ λ₁, λ₂` (what `eig::sym_eig3` returns), then `λ₀` is the minimum of Eq. 3 over unit
  vectors, attained at `e₀` (Eq. 3 = smallest eigenvalue of Eq. 4, Rayleigh);
* `eigenvalue_supergradient` — `λ_min(B) ≤ λ_min(A) + xᵀ(B − A)x` for a unit eigenvector `x` of
  `λ_min(A)`;
* `hasDerivAt_quad_covMatrix` — along differentiable angles `θ_σ(t)` and normals `N_σ(t)`, the
  derivative of `xᵀ A_v(t) x` (fixed `x`) is `Σ (s_σ θ_σ' + y_σ · N_σ')` with Eq. 10's
  `s = (x·N)²`, `y = 2θ(x·N)x`. With the cited derivative of a simple eigenvalue
  (`dev-simple-eigenvalue-derivative`: `dλ = xᵀ dA x`) this is Eq. 10;
* `covMatrix_split_corner` — splitting a corner into two corners with the same normal (angles
  `θ₁ + θ₂ = θ`) leaves `A_v` unchanged.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open scoped DdgProofs.Oracle
open Ddg.Dev (Sym3)

/-! ## Helpers -/

private theorem sumTo_zero' (f : ℕ → ℝ) : Ddg.sumTo 0 f = 0 := by simp [Ddg.sumTo]

private theorem sumTo_succ' (n : ℕ) (f : ℕ → ℝ) : Ddg.sumTo (n + 1) f = Ddg.sumTo n f + f n :=
  rfl

private theorem sumTo_congr (n : ℕ) {f g : ℕ → ℝ} (h : ∀ m < n, f m = g m) :
    Ddg.sumTo n f = Ddg.sumTo n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [sumTo_succ', sumTo_succ', ih (fun m hm => h m (by omega)), h n (by omega)]

/-- Changing one term `j < n` of a sum. -/
private theorem sumTo_update (n j : ℕ) (hj : j < n) {f g : ℕ → ℝ}
    (h : ∀ m < n, m ≠ j → f m = g m) : Ddg.sumTo n f = Ddg.sumTo n g - g j + f j := by
  induction n with
  | zero => omega
  | succ n ih =>
    rw [sumTo_succ', sumTo_succ']
    rcases Nat.lt_succ_iff_lt_or_eq.mp hj with hj' | rfl
    · rw [ih hj' (fun m hm hmj => h m (by omega) hmj), h n (by omega) (by omega)]
      ring
    · rw [sumTo_congr j (fun m hm => h m (by omega) (by omega))]
      ring

theorem sym3_ext {A B : Sym3 ℝ} (h1 : A.xx = B.xx) (h2 : A.yy = B.yy) (h3 : A.zz = B.zz)
    (h4 : A.xy = B.xy) (h5 : A.xz = B.xz) (h6 : A.yz = B.yz) : A = B := by
  cases A; cases B; simp_all

open Ddg.Dev in
theorem quad_add (A B : Sym3 ℝ) (u : V3) : Sym3.quad (Sym3.add A B) u = Sym3.quad A u + Sym3.quad B u := by
  simp only [Sym3.quad, Sym3.apply, Sym3.add, dot_def]; ring

open Ddg.Dev in
theorem quad_zero (u : V3) : Sym3.quad (Sym3.zero : Sym3 ℝ) u = 0 := by
  simp [Sym3.quad, Sym3.apply, Sym3.zero, dot_def]

open Ddg.Dev in
theorem quad_outer (θ : ℝ) (N u : V3) : Sym3.quad (Sym3.outer θ N) u = θ * (dot u N * dot u N) := by
  simp only [Sym3.quad, Sym3.apply, Sym3.outer, dot_def]; ring

section dotlin

open Ddg.Dev

variable (p q r : V3) (c : ℝ)

theorem dot_add_left' : dot (add p q) r = dot p r + dot q r := by simp only [dot_def, add_x, add_y, add_z]; ring
theorem dot_add_right' : dot p (add q r) = dot p q + dot p r := by simp only [dot_def, add_x, add_y, add_z]; ring
theorem dot_smul_left' : dot (smul c p) q = c * dot p q := by simp only [dot_def, smul_x, smul_y, smul_z]; ring
theorem dot_smul_right' : dot p (smul c q) = c * dot p q := by simp only [dot_def, smul_x, smul_y, smul_z]; ring
theorem apply_add (A : Sym3 ℝ) : Sym3.apply A (add p q) = add (Sym3.apply A p) (Sym3.apply A q) := by
  ext <;> simp [Sym3.apply] <;> ring
theorem apply_smul (A : Sym3 ℝ) : Sym3.apply A (smul c p) = smul c (Sym3.apply A p) := by
  ext <;> simp [Sym3.apply] <;> ring

end dotlin

/-! ## Eq. 3 = Eq. 4 -/

open Ddg.Dev in
/-- `uᵀ A_v u = Σ_σ θ_σ ⟨u, N_σ⟩²`: the covariance matrix is the matrix of the objective of Eq. 3. -/
theorem quad_covMatrix (n : ℕ) (θ : ℕ → ℝ) (N : ℕ → V3) (u : V3) :
    Sym3.quad (covMatrix n θ N) u = covObjective n θ N u := by
  induction n with
  | zero => simp [covMatrix, covObjective, Sym3.sumTo, quad_zero, Ddg.sumTo]
  | succ n ih =>
    change Sym3.quad (Sym3.add (covMatrix n θ N) (Sym3.outer (θ n) (N n))) u =
      covObjective n θ N u + θ n * (dot u (N n) * dot u (N n))
    rw [quad_add, ih, quad_outer]

section rayleigh

open Ddg.Dev

variable {A : Sym3 ℝ} {e0 e1 e2 : V3} {l0 l1 l2 : ℝ}

/-- `e₀, e₁, e₂` is an orthonormal basis of `ℝ³` of eigenvectors of `A` with eigenvalues
`l₀, l₁, l₂` (the output of a symmetric eigensolver). -/
structure IsEigenbasis (A : Sym3 ℝ) (e0 e1 e2 : V3) (l0 l1 l2 : ℝ) : Prop where
  n0 : dot e0 e0 = 1
  n1 : dot e1 e1 = 1
  n2 : dot e2 e2 = 1
  o01 : dot e0 e1 = 0
  o02 : dot e0 e2 = 0
  o12 : dot e1 e2 = 0
  a0 : Sym3.apply A e0 = smul l0 e0
  a1 : Sym3.apply A e1 = smul l1 e1
  a2 : Sym3.apply A e2 = smul l2 e2

/-- Completeness of an orthonormal triple of `ℝ³`: `u = Σ ⟨e_i, u⟩ e_i`. -/
theorem expand_orthonormal (h : IsEigenbasis A e0 e1 e2 l0 l1 l2) (u : V3) :
    u = add (add (smul (dot e0 u) e0) (smul (dot e1 u) e1)) (smul (dot e2 u) e2) := by
  let E : Matrix (Fin 3) (Fin 3) ℝ := !![e0.x, e0.y, e0.z; e1.x, e1.y, e1.z; e2.x, e2.y, e2.z]
  have hEE : E * E.transpose = 1 := by
    have n0 := h.n0; have n1 := h.n1; have n2 := h.n2
    have o01 := h.o01; have o02 := h.o02; have o12 := h.o12
    simp only [dot_def] at n0 n1 n2 o01 o02 o12
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [E, Matrix.mul_apply, Fin.sum_univ_three] <;> linarith
  have hEE' : E.transpose * E = 1 := mul_eq_one_comm.mp hEE
  have K : ∀ a b : Fin 3, (E.transpose * E) a b = (1 : Matrix (Fin 3) (Fin 3) ℝ) a b := by
    intro a b; rw [hEE']
  have Kxx := K 0 0; have Kxy := K 0 1; have Kxz := K 0 2
  have Kyy := K 1 1; have Kyz := K 1 2; have Kzz := K 2 2
  simp [E, Matrix.mul_apply, Fin.sum_univ_three] at Kxx Kxy Kxz Kyy Kyz Kzz
  ext <;> simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z, dot_def]
  · linear_combination (-u.x) * Kxx - u.y * Kxy - u.z * Kxz
  · linear_combination (-u.x) * Kxy - u.y * Kyy - u.z * Kyz
  · linear_combination (-u.x) * Kxz - u.y * Kyz - u.z * Kzz

/-- In an orthonormal eigenbasis, `uᵀ A u = Σ λ_i ⟨e_i, u⟩²` and `|u|² = Σ ⟨e_i, u⟩²`. -/
theorem quad_eq_of_eigenbasis (h : IsEigenbasis A e0 e1 e2 l0 l1 l2) (u : V3) :
    Sym3.quad A u = l0 * dot e0 u ^ 2 + l1 * dot e1 u ^ 2 + l2 * dot e2 u ^ 2 ∧
      norm2 u = dot e0 u ^ 2 + dot e1 u ^ 2 + dot e2 u ^ 2 := by
  have hu := expand_orthonormal h u
  set c0 := dot e0 u
  set c1 := dot e1 u
  set c2 := dot e2 u
  have h10 : dot e1 e0 = 0 := by rw [dot_comm]; exact h.o01
  have h20 : dot e2 e0 = 0 := by rw [dot_comm]; exact h.o02
  have h21 : dot e2 e1 = 0 := by rw [dot_comm]; exact h.o12
  constructor
  · rw [Sym3.quad, hu, apply_add, apply_add, apply_smul, apply_smul, apply_smul, h.a0, h.a1, h.a2]
    simp only [dot_add_left', dot_add_right', dot_smul_left', dot_smul_right', h.n0, h.n1, h.n2,
      h.o01, h.o02, h.o12, h10, h20, h21]
    ring
  · rw [norm2, hu]
    simp only [dot_add_left', dot_add_right', dot_smul_left', dot_smul_right', h.n0, h.n1, h.n2,
      h.o01, h.o02, h.o12, h10, h20, h21]
    ring

/-- **Rayleigh**: the smallest eigenvalue is the minimum of the quadratic form over unit vectors,
attained at its eigenvector. -/
theorem quad_isLeast (h : IsEigenbasis A e0 e1 e2 l0 l1 l2) (h01 : l0 ≤ l1) (h02 : l0 ≤ l2) :
    IsLeast {q | ∃ u, norm2 u = 1 ∧ Sym3.quad A u = q} l0 := by
  constructor
  · refine ⟨e0, h.n0, ?_⟩
    have := (quad_eq_of_eigenbasis h e0).1
    rw [this, h.n0, dot_comm e1 e0, dot_comm e2 e0, h.o01, h.o02]; ring
  · rintro q ⟨u, hu, rfl⟩
    obtain ⟨hq, hn⟩ := quad_eq_of_eigenbasis h u
    rw [hq]
    rw [hu] at hn
    have key : l0 * dot e0 u ^ 2 + l1 * dot e1 u ^ 2 + l2 * dot e2 u ^ 2 - l0 =
        (l1 - l0) * dot e1 u ^ 2 + (l2 - l0) * dot e2 u ^ 2 := by
      linear_combination (-l0) * hn
    nlinarith [mul_nonneg (sub_nonneg.mpr h01) (sq_nonneg (dot e1 u)),
      mul_nonneg (sub_nonneg.mpr h02) (sq_nonneg (dot e2 u))]

/-- **Eq. 3 = smallest eigenvalue of Eq. 4** (SPEC §3.4): if `e₀, e₁, e₂` is an orthonormal
eigenbasis of `A_v = Σ θ_σ N_σ N_σᵀ` with eigenvalues `λ₀ ≤ λ₁, λ₂`, then `λ₀` is the minimum of
`Σ_σ θ_σ ⟨u, N_σ⟩²` over unit vectors `u`, attained at `u = e₀`.

[formula: dev-covariance-matrix] -/
theorem covObjective_isLeast (n : ℕ) (θ : ℕ → ℝ) (N : ℕ → V3)
    (h : IsEigenbasis (covMatrix n θ N) e0 e1 e2 l0 l1 l2) (h01 : l0 ≤ l1) (h02 : l0 ≤ l2) :
    IsLeast {q | ∃ u, norm2 u = 1 ∧ covObjective n θ N u = q} l0 := by
  have := quad_isLeast h h01 h02
  simpa only [quad_covMatrix] using this

/-- `A − B`. -/
def sym3Sub (A B : Sym3 ℝ) : Sym3 ℝ :=
  ⟨A.xx - B.xx, A.yy - B.yy, A.zz - B.zz, A.xy - B.xy, A.xz - B.xz, A.yz - B.yz⟩

theorem quad_sub (A B : Sym3 ℝ) (u : V3) :
    Sym3.quad (sym3Sub A B) u = Sym3.quad A u - Sym3.quad B u := by
  simp only [Sym3.quad, Sym3.apply, sym3Sub, dot_def]; ring

/-- **Supergradient of the smallest eigenvalue** (SPEC §3.4): let `x` be a unit vector with
`xᵀ A x = λ_A` (e.g. a unit eigenvector of the smallest eigenvalue `λ_A` of `A`), and let `λ_B`
be a lower bound of `uᵀ B u` on unit vectors (e.g. the smallest eigenvalue of `B`,
`covObjective_isLeast`). Then `λ_B ≤ λ_A + xᵀ(B − A)x`: the concave function `λ_min` lies below
its linearization with `∇_A λ = x xᵀ`, whatever the multiplicity.

[formula: dev-eigenvalue-supergradient] -/
theorem eigenvalue_supergradient (A B : Sym3 ℝ) {x : V3} {lA lB : ℝ} (hx : norm2 x = 1)
    (hA : Sym3.quad A x = lA) (hB : ∀ u, norm2 u = 1 → lB ≤ Sym3.quad B u) :
    lB ≤ lA + Sym3.quad (sym3Sub B A) x := by
  rw [quad_sub, hA]
  have := hB x hx
  linarith

/-- The supergradient inequality between two eigendecompositions: with `λ₀(A)` the smallest
eigenvalue of `A` (eigenvector `e₀`) and `μ₀` that of `B`,
`μ₀ ≤ λ₀ + e₀ᵀ(B − A)e₀`. -/
theorem eigenvalue_supergradient_of_eigenbasis {A B : Sym3 ℝ} {f0 f1 f2 : V3} {m0 m1 m2 : ℝ}
    (hA : IsEigenbasis A e0 e1 e2 l0 l1 l2) (hB : IsEigenbasis B f0 f1 f2 m0 m1 m2)
    (hm1 : m0 ≤ m1) (hm2 : m0 ≤ m2) :
    m0 ≤ l0 + Sym3.quad (sym3Sub B A) e0 := by
  refine eigenvalue_supergradient A B hA.n0 ?_ ?_
  · have := (quad_eq_of_eigenbasis hA e0).1
    rw [this, hA.n0, dot_comm e1 e0, dot_comm e2 e0, hA.o01, hA.o02]; ring
  · intro u hu
    exact (quad_isLeast hB hm1 hm2).2 ⟨u, hu, rfl⟩

end rayleigh

/-! ## Eq. 10 -/

open Ddg.Dev in
/-- **Gradient of the covariance energy** (Eq. 10, SPEC §3.4): along differentiable corner angles
`θ_σ(t)` and normals `N_σ(t)`, for a fixed vector `x`,
`d/dt xᵀ A_v(t) x = Σ_σ (s_σ θ_σ' + y_σ · N_σ')` with `s_σ = (x·N_σ)²` (`covSensitivityTheta`)
and `y_σ = 2θ_σ (x·N_σ) x` (`covSensitivityN`). For `x` the unit eigenvector of a simple smallest
eigenvalue, `dλ = xᵀ dA x` (cited, `dev-simple-eigenvalue-derivative`), so this is `dλ`; composed
with the face derivatives of Eq. 8, 9 (`dev-normal-gradient`, `dev-angle-gradient`) it is Eq. 10.

[formula: dev-covariance-gradient] -/
theorem hasDerivAt_quad_covMatrix (n : ℕ) (x : V3) (θ : ℕ → ℝ → ℝ) (θ' : ℕ → ℝ)
    (N : ℕ → ℝ → V3) (N' : ℕ → V3) (t : ℝ)
    (hθ : ∀ m < n, HasDerivAt (θ m) (θ' m) t) (hN : ∀ m < n, HasDerivAt3 (N m) (N' m) t) :
    HasDerivAt (fun s => Sym3.quad (covMatrix n (fun m => θ m s) (fun m => N m s)) x)
      (Ddg.sumTo n (fun m => covSensitivityTheta x (N m t) * θ' m +
        dot (covSensitivityN (θ m t) x (N m t)) (N' m))) t := by
  simp only [quad_covMatrix]
  induction n with
  | zero =>
    simp only [covObjective, sumTo_zero']
    exact hasDerivAt_const t 0
  | succ n ih =>
    have ih' := ih (fun m hm => hθ m (by omega)) (fun m hm => hN m (by omega))
    simp only [covObjective, sumTo_succ'] at ih' ⊢
    have hd := (hN n (by omega)).dot_const x
    have hterm := (hθ n (by omega)).mul (hd.mul hd)
    convert ih'.add hterm using 1
    simp only [covSensitivityTheta, covSensitivityN, dot_smul_left', Pi.mul_apply]
    ring

/-! ## Retessellation -/

open Ddg.Dev in
/-- **Angle weights and retessellation** (SPEC §3.4): splitting the corner `j` of a star into two
corners with the same normal — the new star has the corners `0 … n` with `θ'_j + θ'_n = θ_j`,
`N'_n = N_j`, and all other corners unchanged — leaves `A_v = Σ θ N Nᵀ` unchanged. (Splitting an
edge of a planar region at a point does exactly this to the corner it crosses.)

[formula: dev-angle-weight-retessellation] -/
theorem covMatrix_split_corner (n j : ℕ) (hj : j < n) (θ θ' : ℕ → ℝ) (N N' : ℕ → V3)
    (hθ : ∀ m < n, m ≠ j → θ' m = θ m) (hN : ∀ m < n, N' m = N m) (hNn : N' n = N j)
    (hsplit : θ' j + θ' n = θ j) :
    covMatrix (n + 1) θ' N' = covMatrix n θ N := by
  -- every entry is a sum `Σ θ_m E(N_m)`
  have entry : ∀ (E : V3 → ℝ) (P : Sym3 ℝ → ℝ),
      (∀ k (θ'' : ℕ → ℝ) (M : ℕ → V3), P (covMatrix k θ'' M) =
        Ddg.sumTo k (fun m => θ'' m * E (M m))) →
      P (covMatrix (n + 1) θ' N') = P (covMatrix n θ N) := by
    intro E P hP
    rw [hP, hP, sumTo_succ',
      sumTo_update n j hj (g := fun m => θ m * E (N m)) (fun m hm hmj => by
        simp only [hθ m hm hmj, hN m hm])]
    simp only [hN j hj, hNn]
    linear_combination E (N j) * hsplit
  have sumEntry : ∀ (P : Sym3 ℝ → ℝ) (hadd : ∀ A B, P (Sym3.add A B) = P A + P B)
      (hzero : P Sym3.zero = 0) (E : V3 → ℝ) (hout : ∀ θ M, P (Sym3.outer θ M) = θ * E M),
      ∀ k (θ'' : ℕ → ℝ) (M : ℕ → V3), P (covMatrix k θ'' M) =
        Ddg.sumTo k (fun m => θ'' m * E (M m)) := by
    intro P hadd hzero E hout k θ'' M
    induction k with
    | zero => simp [covMatrix, Sym3.sumTo, hzero, sumTo_zero']
    | succ k ih =>
      change P (Sym3.add (covMatrix k θ'' M) (Sym3.outer (θ'' k) (M k))) = _
      rw [hadd, ih, hout, sumTo_succ']
  apply sym3_ext
  · exact entry (fun M => M.x * M.x) Sym3.xx
      (sumEntry Sym3.xx (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))
  · exact entry (fun M => M.y * M.y) Sym3.yy
      (sumEntry Sym3.yy (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))
  · exact entry (fun M => M.z * M.z) Sym3.zz
      (sumEntry Sym3.zz (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))
  · exact entry (fun M => M.x * M.y) Sym3.xy
      (sumEntry Sym3.xy (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))
  · exact entry (fun M => M.x * M.z) Sym3.xz
      (sumEntry Sym3.xz (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))
  · exact entry (fun M => M.y * M.z) Sym3.yz
      (sumEntry Sym3.yz (fun _ _ => rfl) (by simp [Sym3.zero]) _ (fun _ _ => rfl))

end Dev

end DdgProofs
