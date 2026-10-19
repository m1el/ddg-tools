import DdgProofs.Developability.Vec

/-!
# Developability: the intrinsic width — log map, its adjoint, the tangent eigenvalue (SPEC §3.6)

`crates/ddg-developability/SPEC.md` §3.6, decisions D6, D12 (paper App. B.5, Eq. 11).
With the vertex normal `n0` and a face normal `N` (both unit), `c = ⟨n0, N⟩`, `s = |n0 × N|`,
`φ = atan2(s, c)`:

* `logMap_spec` (`dev-log-map`): `cos φ = c`, `sin φ = s`, `Ñ = log_{n0}(N) = (φ/s)(N − c n0)` is
  tangent at `n0` (`Ñ ⊥ n0`), has length `φ`, and `exp_{n0}(Ñ) = N`; `logMap_of_parallel`: for
  `s = 0`, `c > 0` the normals agree and `Ñ = 0`.
* `hasDerivAt_logMap` (`dev-log-map-gradient`): along curves `n0(t)`, `N(t)` that are unit at `t`
  with velocities tangent to the sphere (`n0·n0' = 0`, `N·N' = 0`, as for curves of unit vectors),
  `d/dt ⟨z, Ñ⟩ = y_N·N' + y_{n0}·n0'` with `(y_N, y_{n0}) = logMapAdjoint n0 N z`:
  `K = G (z·n0) + H z·(N − c n0)`, `y_N = G z − K n0`, `y_{n0} = −G c z − K N`, `G = φ/s`,
  `H = (s − φ c)/s³`. This is the paper's Eq. 11 with the correction D6 (`μ_i μ_iᵀ` in the second
  bracket), written in the smooth form of the SPEC.
* `hasDerivAt_quad_tangent` (`dev-tangent-eigen-derivative`): the derivative of
  `yᵀ T(t)ᵀ Ã(t) T(t) y` (the 2x2 Rayleigh quotient in a moving orthonormal tangent basis
  `T = [t1 t2]`, fixed `y`) is `xᵀ Ã' x`, `x = T y`, whenever `Ã x = λ x`: the terms of the moving
  basis vanish. `apply_eq_smul_of_tangent_eigen`: an eigenvector `y` of the 2x2 matrix `TᵀÃT` gives
  the eigenvector `T y` of the 3x3 matrix `Ã = Σ θ Ñ Ñᵀ` when all `Ñ` lie in the plane of `t1, t2`.
  With the derivative of a simple eigenvalue (`dev-simple-eigenvalue-derivative`, cited) this is
  `dλ = xᵀ dÃ x` with the 3D `dÃ`.

The exact coefficients `G`, `H` are covered; their Taylor polynomials near `φ = 0`
(`logGSeries`, `logHSeries`) are a numerical device and are not.
-/

noncomputable section

namespace DdgProofs

namespace Dev

open scoped DdgProofs.Oracle
open Real Ddg.Dev

/-! ## The log map -/

section logmap

variable (n0 N : V3)

/-- `s² + c² = 1` for unit `n0`, `N`. -/
theorem norm_cross_sq_add_dot_sq (hn0 : norm2 n0 = 1) (hN : norm2 N = 1) :
    dot n0 N ^ 2 + norm (cross n0 N) ^ 2 = 1 := by
  have := norm2_cross_add_sq_dot n0 N
  rw [hn0, hN, ← norm_sq] at this
  linarith

/-- **The log map** (SPEC §3.6, paper App. B.5): for unit `n0`, `N` with `s = |n0 × N| > 0`,
`c = ⟨n0, N⟩` and `φ = atan2(s, c)`: `cos φ = c`, `sin φ = s`, `0 < φ`, the vector
`Ñ = logMap n0 N = (φ/s)(N − c n0)` is orthogonal to `n0`, has length `φ`, and
`exp_{n0}(Ñ) = cos φ · n0 + (sin φ/φ) Ñ = N`.

[formula: dev-log-map] -/
theorem logMap_spec (hn0 : norm2 n0 = 1) (hN : norm2 N = 1) (hs : 0 < norm (cross n0 N)) :
    cos (logAngle n0 N) = dot n0 N ∧ sin (logAngle n0 N) = norm (cross n0 N) ∧
      0 < logAngle n0 N ∧ dot n0 (logMap n0 N) = 0 ∧ norm (logMap n0 N) = logAngle n0 N ∧
      expMapWith (cos (logAngle n0 N)) (sin (logAngle n0 N) / logAngle n0 N) n0 (logMap n0 N) =
        N := by
  have hcs := norm_cross_sq_add_dot_sq n0 N hn0 hN
  have hpos : 0 < dot n0 N ^ 2 + norm (cross n0 N) ^ 2 := by rw [hcs]; norm_num
  have hcos := cos_atan2 (norm_nonneg' (cross n0 N)) hpos
  have hsin := sin_atan2 (norm_nonneg' (cross n0 N)) hpos
  rw [hcs, sqrt_one, div_one] at hcos hsin
  change cos (logAngle n0 N) = _ at hcos
  change sin (logAngle n0 N) = _ at hsin
  have hφ : 0 < logAngle n0 N := by
    show 0 < Ddg.Dev.atan2 _ _
    rw [atan2_of_pos hs]
    linarith [arctan_lt_pi_div_two (dot n0 N / norm (cross n0 N))]
  set s := norm (cross n0 N) with hsdef
  set c := dot n0 N with hcdef
  set φ := logAngle n0 N with hφdef
  have hlog : logMap n0 N = smul (φ / s) (sub N (smul c n0)) := rfl
  have hperp : dot n0 (sub N (smul c n0)) = 0 := by
    have : dot n0 (sub N (smul c n0)) = c - c * norm2 n0 := by
      simp only [hcdef, dot_def, norm2_def, sub_x, sub_y, sub_z, smul_x, smul_y, smul_z]; ring
    rw [this, hn0]; ring
  have hv2 : norm2 (sub N (smul c n0)) = s ^ 2 := by
    have : norm2 (sub N (smul c n0)) = norm2 N - 2 * c * c + c ^ 2 * norm2 n0 := by
      simp only [hcdef, dot_def, norm2_def, sub_x, sub_y, sub_z, smul_x, smul_y, smul_z]; ring
    rw [this, hn0, hN]; linarith
  have hvn : norm (sub N (smul c n0)) = s := by
    rw [norm_def, hv2, sqrt_sq (norm_nonneg' _)]
  refine ⟨hcos, hsin, hφ, ?_, ?_, ?_⟩
  · rw [hlog]
    have : dot n0 (smul (φ / s) (sub N (smul c n0))) = φ / s * dot n0 (sub N (smul c n0)) := by
      simp only [dot_def, smul_x, smul_y, smul_z]; ring
    rw [this, hperp, mul_zero]
  · rw [hlog, norm_smul_of_nonneg _ (div_nonneg hφ.le hs.le), hvn]
    field_simp
  · rw [hcos, hsin, hlog]
    unfold expMapWith
    ext <;> simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z, sub_x, sub_y, sub_z] <;>
      field_simp <;> ring

/-- **The log map at `φ = 0`**: for unit `n0`, `N` with `n0 × N = 0` and `⟨n0, N⟩ > 0`, the
normals agree and `logMap n0 N = 0` (`exp_{n0}(0) = n0 = N`). -/
theorem logMap_of_parallel (hn0 : norm2 n0 = 1) (hN : norm2 N = 1)
    (hs : norm (cross n0 N) = 0) (hc : 0 < dot n0 N) : N = n0 ∧ logMap n0 N = zero := by
  have hcs := norm_cross_sq_add_dot_sq n0 N hn0 hN
  rw [hs] at hcs
  have hc1 : dot n0 N = 1 := by nlinarith
  have hd : norm2 (sub N n0) = 0 := by
    have : norm2 (sub N n0) = norm2 N - 2 * dot n0 N + norm2 n0 := by
      simp only [dot_def, norm2_def, sub_x, sub_y, sub_z]; ring
    rw [this, hN, hn0, hc1]; ring
  rw [norm2_eq_zero_iff] at hd
  have hNn : N = n0 := by
    ext
    · have := congrArg Ddg.Point.x hd; simp at this; linarith
    · have := congrArg Ddg.Point.y hd; simp at this; linarith
    · have := congrArg Ddg.Point.z hd; simp at this; linarith
  refine ⟨hNn, ?_⟩
  have : logMap n0 N = smul (logG n0 N) (sub N (smul (dot n0 N) n0)) := rfl
  rw [this, hc1, hNn]
  ext <;> simp

end logmap

/-! ## The adjoint of the log map -/

/-- Binet–Cauchy for the derivative of `|a × b|²`. -/
theorem dot_cross_deriv (a b a' b' : V3) :
    dot (cross a b) (add (cross a' b) (cross a b')) =
      dot a a' * norm2 b + norm2 a * dot b b' - dot a b * (dot a' b + dot a b') := by
  simp only [dot_def, norm2_def, add_x, add_y, add_z, cross_x, cross_y, cross_z]; ring

/-- **Derivative of the log map; its adjoint** (SPEC §3.6, corrected Eq. 11, D6): let `n0(t)`,
`N(t)` be differentiable curves with velocities `n0'`, `N'` that are unit at `t` with tangent
velocities (`n0·n0' = 0`, `N·N' = 0`; true for curves of unit vectors) and `|n0 × N| > 0` at `t`
(`0 < φ < π`). Then for every sensitivity `z`,
`d/dt ⟨z, log_{n0}(N)⟩ = y_N·N' + y_{n0}·n0'` with `(y_N, y_{n0}) = logMapAdjoint n0 N z`
(`K = G (z·n0) + H z·(N − c n0)`, `y_N = G z − K n0`, `y_{n0} = −G c z − K N`).

[formula: dev-log-map-gradient] -/
theorem hasDerivAt_logMap {n0 N : ℝ → V3} {n0' N' : V3} {t : ℝ} (z : V3)
    (h1 : HasDerivAt3 n0 n0' t) (h2 : HasDerivAt3 N N' t) (hn0 : norm2 (n0 t) = 1)
    (hN : norm2 (N t) = 1) (hn0' : dot (n0 t) n0' = 0) (hN' : dot (N t) N' = 0)
    (hs : 0 < norm (cross (n0 t) (N t))) :
    HasDerivAt (fun s => dot z (logMap (n0 s) (N s)))
      (dot (logMapAdjoint (n0 t) (N t) z).1 N' + dot (logMapAdjoint (n0 t) (N t) z).2 n0') t := by
  have hX := h1.cross h2
  have hX0 : cross (n0 t) (N t) ≠ zero := (norm_pos_iff' _).mp hs
  have hc := h1.dot h2
  set c' := dot n0' (N t) + dot (n0 t) N' with hc'
  have hcs := norm_cross_sq_add_dot_sq (n0 t) (N t) hn0 hN
  -- `s' = −c c'/s`
  have hsn : HasDerivAt (fun s => norm (cross (n0 s) (N s)))
      (-(dot (n0 t) (N t) * c') / norm (cross (n0 t) (N t))) t := by
    refine (hX.norm hX0).congr_deriv ?_
    rw [dot_cross_deriv, hn0, hN, hn0', hN']
    ring
  -- `φ' = −c'/s`
  have hφ : HasDerivAt (fun s => logAngle (n0 s) (N s)) (-c' / norm (cross (n0 t) (N t))) t := by
    refine (Dev.HasDerivAt.atan2 hc hsn hs).congr_deriv ?_
    rw [hcs]
    field_simp
    linear_combination (-c') * hcs
  -- `G' = −H c'`
  have hG : HasDerivAt (fun s => logG (n0 s) (N s)) (-(logH (n0 t) (N t) * c')) t := by
    refine (hφ.div hsn hs.ne').congr_deriv ?_
    simp only [logH]
    field_simp
    ring
  have hV := h2.sub (HasDerivAt3.smul hc h1)
  have hL := (HasDerivAt3.smul hG hV).dot_const z
  refine hL.congr_deriv ?_
  simp only [logMapAdjoint, logMapAdjointWith, hc', dot_def, add_x, add_y, add_z, sub_x, sub_y,
    sub_z, smul_x, smul_y, smul_z, neg_x, neg_y, neg_z]
  ring

/-! ## The tangent eigenvalue -/

/-- `t` has the derivative `A'` at `t`, entry by entry. -/
structure HasDerivAtSym (A : ℝ → Sym3 ℝ) (A' : Sym3 ℝ) (t : ℝ) : Prop where
  xx : HasDerivAt (fun s => (A s).xx) A'.xx t
  yy : HasDerivAt (fun s => (A s).yy) A'.yy t
  zz : HasDerivAt (fun s => (A s).zz) A'.zz t
  xy : HasDerivAt (fun s => (A s).xy) A'.xy t
  xz : HasDerivAt (fun s => (A s).xz) A'.xz t
  yz : HasDerivAt (fun s => (A s).yz) A'.yz t

theorem Sym3.apply_x (A : Sym3 ℝ) (v : V3) :
    (Sym3.apply A v).x = A.xx * v.x + A.xy * v.y + A.xz * v.z := rfl
theorem Sym3.apply_y (A : Sym3 ℝ) (v : V3) :
    (Sym3.apply A v).y = A.xy * v.x + A.yy * v.y + A.yz * v.z := rfl
theorem Sym3.apply_z (A : Sym3 ℝ) (v : V3) :
    (Sym3.apply A v).z = A.xz * v.x + A.yz * v.y + A.zz * v.z := rfl

/-- Derivative of `A(t) v(t)`. -/
theorem HasDerivAt3.sym3_apply {A : ℝ → Sym3 ℝ} {A' : Sym3 ℝ} {v : ℝ → V3} {v' : V3} {t : ℝ}
    (hA : HasDerivAtSym A A' t) (hv : HasDerivAt3 v v' t) :
    HasDerivAt3 (fun s => Sym3.apply (A s) (v s))
      (Ddg.Dev.add (Sym3.apply A' (v t)) (Sym3.apply (A t) v')) t := by
  refine ⟨?_, ?_, ?_⟩
  · convert ((hA.xx.mul hv.x).add (hA.xy.mul hv.y)).add (hA.xz.mul hv.z) using 1 <;>
      first | rfl | (simp only [add_x, Sym3.apply_x]; ring)
  · convert ((hA.xy.mul hv.x).add (hA.yy.mul hv.y)).add (hA.yz.mul hv.z) using 1 <;>
      first | rfl | (simp only [add_y, Sym3.apply_y]; ring)
  · convert ((hA.xz.mul hv.x).add (hA.yz.mul hv.y)).add (hA.zz.mul hv.z) using 1 <;>
      first | rfl | (simp only [add_z, Sym3.apply_z]; ring)

/-- `yᵀ (TᵀÃT) y = xᵀ Ã x` with `x = T y`: the 2x2 Rayleigh quotient is the 3x3 one at `T y`. -/
theorem quad_tangent (A : Sym3 ℝ) (t1 t2 : V3) (y1 y2 : ℝ) :
    y1 * y1 * dot t1 (Sym3.apply A t1) + 2 * (y1 * y2) * dot t1 (Sym3.apply A t2) +
        y2 * y2 * dot t2 (Sym3.apply A t2) =
      Sym3.quad A (add (smul y1 t1) (smul y2 t2)) := by
  simp only [Sym3.quad, dot_def, Sym3.apply_x, Sym3.apply_y, Sym3.apply_z, add_x, add_y, add_z,
    smul_x, smul_y, smul_z]
  ring

/-- **The moving tangent basis drops out of the eigenvalue derivative** (SPEC §3.6): let
`t1(t)`, `t2(t)` be differentiable curves, orthonormal to first order at `t`
(`t1·t1' = 0`, `t2·t2' = 0`, `t1'·t2 + t1·t2' = 0`, as for curves of orthonormal pairs), `Ã(t)` a
differentiable curve of symmetric matrices, `y = (y1, y2)` fixed, `x(t) = y1 t1(t) + y2 t2(t)`, and
`Ã x = λ x` at `t`. Then `d/dt [yᵀ T(t)ᵀ Ã(t) T(t) y] = xᵀ Ã' x`: the terms
`yᵀ(dTᵀ Ã T + Tᵀ Ã dT) y = 2λ xᵀ dT y` vanish.

[formula: dev-tangent-eigen-derivative] -/
theorem hasDerivAt_quad_tangent {A : ℝ → Sym3 ℝ} {A' : Sym3 ℝ} {t1 t2 : ℝ → V3} {t1' t2' : V3}
    {t : ℝ} (y1 y2 lam : ℝ) (hA : HasDerivAtSym A A' t) (h1 : HasDerivAt3 t1 t1' t)
    (h2 : HasDerivAt3 t2 t2' t) (o11 : dot (t1 t) t1' = 0) (o22 : dot (t2 t) t2' = 0)
    (o12 : dot t1' (t2 t) + dot (t1 t) t2' = 0)
    (heig : Sym3.apply (A t) (add (smul y1 (t1 t)) (smul y2 (t2 t))) =
      smul lam (add (smul y1 (t1 t)) (smul y2 (t2 t)))) :
    HasDerivAt (fun s => Sym3.quad (A s) (add (smul y1 (t1 s)) (smul y2 (t2 s))))
      (Sym3.quad A' (add (smul y1 (t1 t)) (smul y2 (t2 t)))) t := by
  have hx := (h1.smul_const y1).add (h2.smul_const y2)
  have hd := hx.dot (HasDerivAt3.sym3_apply hA hx)
  refine hd.congr_deriv ?_
  set x := add (smul y1 (t1 t)) (smul y2 (t2 t)) with hxdef
  set x' := add (smul y1 t1') (smul y2 t2')
  -- `x'ᵀ Ã x + xᵀ Ã x' = 2 λ x·x' = 0`
  have hsym : dot x (Sym3.apply (A t) x') = dot x' (Sym3.apply (A t) x) := by
    simp only [dot_def, Sym3.apply_x, Sym3.apply_y, Sym3.apply_z]; ring
  have hxx : dot x' x = 0 := by
    have : dot x' x = y1 * y1 * dot (t1 t) t1' + y2 * y2 * dot (t2 t) t2' +
        y1 * y2 * (dot t1' (t2 t) + dot (t1 t) t2') := by
      simp only [x, x', dot_def, add_x, add_y, add_z, smul_x, smul_y, smul_z]; ring
    rw [this, o11, o22, o12]; ring
  have e : dot x (add (Sym3.apply A' x) (Sym3.apply (A t) x')) =
      Sym3.quad A' x + dot x (Sym3.apply (A t) x') := by
    simp only [Sym3.quad, dot_def, add_x, add_y, add_z]; ring
  rw [e, hsym, heig]
  have : dot x' (smul lam x) = lam * dot x' x := by
    simp only [dot_def, smul_x, smul_y, smul_z]; ring
  rw [this, hxx]
  ring

/-! ### The 3D eigenvector from the 2x2 one -/

/-- `(a × b) × a = |a|² b − (a·b) a`. -/
theorem cross_cross_left (a b : V3) :
    cross (cross a b) a = sub (smul (norm2 a) b) (smul (dot a b) a) := by
  ext <;> simp [dot_def, norm2_def] <;> ring

/-- Completeness of the orthonormal frame `(t1, t2, t1 × t2)`. -/
theorem frame_expand_tangent (t1 t2 X : V3) (h1 : norm2 t1 = 1) (h2 : norm2 t2 = 1)
    (h12 : dot t1 t2 = 0) :
    X = add (add (smul (dot t1 X) t1) (smul (dot t2 X) t2))
      (smul (dot (cross t1 t2) X) (cross t1 t2)) := by
  rw [norm2_def] at h1 h2
  rw [dot_def] at h12
  -- Cramer: `det · X = (X·(b×c)) a + (X·(c×a)) b + (X·(a×b)) c` for `a = t1, b = t2, c = t1×t2`
  ext <;> simp only [add_x, add_y, add_z, smul_x, smul_y, smul_z, cross_x, cross_y, cross_z,
      dot_def]
  · linear_combination (exp := 1)
      (-X.x + (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t2.x) * h1 +
      (-(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z) * X.x +
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t1.x) * h2 +
      ((t1.x * t2.x + t1.y * t2.y + t1.z * t2.z) * X.x -
        (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t1.x -
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t2.x) * h12
  · linear_combination (exp := 1)
      (-X.y + (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t2.y) * h1 +
      (-(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z) * X.y +
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t1.y) * h2 +
      ((t1.x * t2.x + t1.y * t2.y + t1.z * t2.z) * X.y -
        (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t1.y -
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t2.y) * h12
  · linear_combination (exp := 1)
      (-X.z + (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t2.z) * h1 +
      (-(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z) * X.z +
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t1.z) * h2 +
      ((t1.x * t2.x + t1.y * t2.y + t1.z * t2.z) * X.z -
        (t2.x * X.x + t2.y * X.y + t2.z * X.z) * t1.z -
        (t1.x * X.x + t1.y * X.y + t1.z * X.z) * t2.z) * h12

/-- `w · (Σ θ Ñ Ñᵀ) v = Σ θ (Ñ·v)(Ñ·w)`. -/
theorem dot_apply_covMatrix (w v : V3) (θ : ℕ → ℝ) (Nt : ℕ → V3) :
    ∀ n, dot w (Sym3.apply (covMatrix n θ Nt) v) =
      Ddg.sumTo n (fun m => θ m * (dot (Nt m) v * dot (Nt m) w))
  | 0 => by
    simp [covMatrix, Ddg.Dev.Sym3.sumTo, Ddg.Dev.Sym3.zero, Ddg.sumTo, dot_def, Sym3.apply_x,
      Sym3.apply_y, Sym3.apply_z]
  | n + 1 => by
    have ih := dot_apply_covMatrix w v θ Nt n
    have e : covMatrix (n + 1) θ Nt = Sym3.add (covMatrix n θ Nt) (Sym3.outer (θ n) (Nt n)) := rfl
    have e2 : Ddg.sumTo (n + 1) (fun m => θ m * (dot (Nt m) v * dot (Nt m) w)) =
        Ddg.sumTo n (fun m => θ m * (dot (Nt m) v * dot (Nt m) w)) +
          θ n * (dot (Nt n) v * dot (Nt n) w) := rfl
    rw [e, e2, ← ih]
    simp only [Sym3.add, Sym3.outer, dot_def, Sym3.apply_x, Sym3.apply_y, Sym3.apply_z]
    ring

theorem sumTo_zero_of (n : ℕ) (f : ℕ → ℝ) (h : ∀ m, f m = 0) : Ddg.sumTo n f = 0 := by
  induction n with
  | zero => simp [Ddg.sumTo]
  | succ n ih => show Ddg.sumTo n f + f n = 0; rw [ih, h]; ring

/-- **The 3D eigenvector from the 2x2 one** (SPEC §3.6): let `t1, t2` be orthonormal, every `Ñ_m`
in their plane (`Ñ_m ⊥ t1 × t2`; in the code `t1 × t2 = ±n0` and `Ñ ⊥ n0`), and `y` an
eigenvector of the 2x2 matrix `TᵀÃT` with eigenvalue `λ` (`Ã = Σ θ Ñ Ñᵀ`). Then `x = y1 t1 + y2 t2`
satisfies `Ã x = λ x` in 3D. -/
theorem apply_eq_smul_of_tangent_eigen (n : ℕ) (θ : ℕ → ℝ) (Nt : ℕ → V3) (t1 t2 : V3)
    (y1 y2 lam : ℝ) (h1 : norm2 t1 = 1) (h2 : norm2 t2 = 1) (h12 : dot t1 t2 = 0)
    (hplane : ∀ m, dot (Nt m) (cross t1 t2) = 0)
    (e1 : dot t1 (Sym3.apply (covMatrix n θ Nt) t1) * y1 +
      dot t1 (Sym3.apply (covMatrix n θ Nt) t2) * y2 = lam * y1)
    (e2 : dot t2 (Sym3.apply (covMatrix n θ Nt) t1) * y1 +
      dot t2 (Sym3.apply (covMatrix n θ Nt) t2) * y2 = lam * y2) :
    Sym3.apply (covMatrix n θ Nt) (add (smul y1 t1) (smul y2 t2)) =
      smul lam (add (smul y1 t1) (smul y2 t2)) := by
  set A := covMatrix n θ Nt
  set x := add (smul y1 t1) (smul y2 t2)
  have hexp := frame_expand_tangent t1 t2 (Sym3.apply A x) h1 h2 h12
  have hn : dot (cross t1 t2) (Sym3.apply A x) = 0 := by
    rw [dot_apply_covMatrix]
    exact sumTo_zero_of _ _ fun m => by rw [hplane m]; ring
  have lin : ∀ w, dot w (Sym3.apply A x) =
      y1 * dot w (Sym3.apply A t1) + y2 * dot w (Sym3.apply A t2) := by
    intro w
    simp only [x, dot_def, Sym3.apply_x, Sym3.apply_y, Sym3.apply_z, add_x, add_y, add_z,
      smul_x, smul_y, smul_z]
    ring
  rw [hn, lin t1, lin t2] at hexp
  rw [hexp]
  have c1 : y1 * dot t1 (Sym3.apply A t1) + y2 * dot t1 (Sym3.apply A t2) = lam * y1 := by
    linarith
  have c2 : y1 * dot t2 (Sym3.apply A t1) + y2 * dot t2 (Sym3.apply A t2) = lam * y2 := by
    linarith
  rw [c1, c2]
  ext <;> simp [x] <;> ring

end Dev

end DdgProofs
