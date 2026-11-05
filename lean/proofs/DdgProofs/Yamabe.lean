import DdgProofs.GaussBonnet

/-!
# The linearized discrete Yamabe equation

Under a discrete conformal change of metric `ℓ̃_ij = e^{(u_i + u_j)/2} ℓ_ij` the corner angles
change, to first order in `u`, by the cotan weights: for a triangle `ijk`,

`dθ_i = ½ cot θ_k (du_j − du_i) + ½ cot θ_j (du_k − du_i)`

(Springborn, Schröder, Pinkall 2008, "Conformal equivalence of triangle meshes", Sec. 4.1/App. B;
Sawhney, Crane 2017, "Boundary First Flattening", Sec. 4.3). Summed over the faces at a vertex:
the derivative of the angle sum with respect to `u` is `−L`, so the discrete curvature
`2π − Σθ` of the scaled metric is `Ω + L u + o(u)`. Prescribing the curvature `μ` gives the
constraint `L u = μ − Ω` of the canonical problem, to first order.
-/

noncomputable section

namespace DDG

open Real

/-- **First variation of a corner angle under a conformal change of the edge lengths.**
Let `a b c` be the side lengths of a nondegenerate triangle (opposite the vertices `i j k`), and
scale them by `e^{t(y+z)/2}`, `e^{t(z+x)/2}`, `e^{t(x+y)/2}` (log scale factors `t·x`, `t·y`, `t·z`
at `i`, `j`, `k`). Then the angle at `i` has derivative
`½ cot θ_k·(y − x) + ½ cot θ_j·(z − x)` at `t = 0`.

[formula: discrete-yamabe-linearization] -/
theorem hasDerivAt_cornerAngle_conformal {a b c : ℝ} (h : IsTriangle a b c) (x y z : ℝ) :
    HasDerivAt
      (fun t : ℝ => cornerAngle (a * exp (t * ((y + z) / 2))) (b * exp (t * ((z + x) / 2)))
        (c * exp (t * ((x + y) / 2))))
      (1 / 2 * cornerCot c a b * (y - x) + 1 / 2 * cornerCot b c a * (z - x)) 0 := by
  have ha := h.pos_a
  have hb := h.pos_b
  have hc := h.pos_c
  have hT := triArea_pos h
  have hbd := cos_arg_bounds h
  set α := (y + z) / 2 with hα
  set β := (z + x) / 2 with hβ
  set γ := (x + y) / 2 with hγ
  have dA : HasDerivAt (fun t : ℝ => a * exp (t * α)) (a * α) 0 := by
    have := ((hasDerivAt_mul_const α (x := (0 : ℝ))).exp).const_mul a
    simpa using this
  have dB : HasDerivAt (fun t : ℝ => b * exp (t * β)) (b * β) 0 := by
    have := ((hasDerivAt_mul_const β (x := (0 : ℝ))).exp).const_mul b
    simpa using this
  have dC : HasDerivAt (fun t : ℝ => c * exp (t * γ)) (c * γ) 0 := by
    have := ((hasDerivAt_mul_const γ (x := (0 : ℝ))).exp).const_mul c
    simpa using this
  have dN : HasDerivAt
      (fun t : ℝ => (b * exp (t * β)) ^ 2 + (c * exp (t * γ)) ^ 2 - (a * exp (t * α)) ^ 2)
      (2 * b * (b * β) + 2 * c * (c * γ) - 2 * a * (a * α)) 0 := by
    have := ((dB.fun_pow 2).fun_add (dC.fun_pow 2)).fun_sub (dA.fun_pow 2)
    convert this using 1
    simp
  have dD : HasDerivAt (fun t : ℝ => 2 * (b * exp (t * β)) * (c * exp (t * γ)))
      (2 * (b * β) * c + 2 * b * (c * γ)) 0 := by
    have := (dB.const_mul 2).fun_mul dC
    convert this using 1
    simp
  have hD0 : (2 * (b * exp (0 * β)) * (c * exp (0 * γ))) ≠ 0 := by
    simp only [zero_mul, exp_zero, mul_one]
    positivity
  have dX : HasDerivAt
      (fun t : ℝ => ((b * exp (t * β)) ^ 2 + (c * exp (t * γ)) ^ 2 - (a * exp (t * α)) ^ 2) /
        (2 * (b * exp (t * β)) * (c * exp (t * γ))))
      (((2 * b * (b * β) + 2 * c * (c * γ) - 2 * a * (a * α)) * (2 * b * c) -
        (b ^ 2 + c ^ 2 - a ^ 2) * (2 * (b * β) * c + 2 * b * (c * γ))) / (2 * b * c) ^ 2) 0 := by
    have := dN.fun_div dD hD0
    convert this using 1
    simp
  have hX0 : ((b * exp (0 * β)) ^ 2 + (c * exp (0 * γ)) ^ 2 - (a * exp (0 * α)) ^ 2) /
      (2 * (b * exp (0 * β)) * (c * exp (0 * γ))) = (b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c) := by
    simp
  have harc : HasDerivAt arccos (-(1 / √(1 - ((b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c)) ^ 2)))
      (((b * exp (0 * β)) ^ 2 + (c * exp (0 * γ)) ^ 2 - (a * exp (0 * α)) ^ 2) /
        (2 * (b * exp (0 * β)) * (c * exp (0 * γ)))) := by
    rw [hX0]
    exact hasDerivAt_arccos hbd.1.ne' hbd.2.ne
  have hcomp := HasDerivAt.comp (0 : ℝ) harc dX
  have hs : √(1 - ((b ^ 2 + c ^ 2 - a ^ 2) / (2 * b * c)) ^ 2) = 2 * triArea a b c / (b * c) := by
    rw [← sin_arccos]
    exact sin_cornerAngle h
  refine hcomp.congr_deriv ?_
  rw [hs]
  unfold cornerCot
  rw [triArea_rotate a b c, triArea_rotate b c a, triArea_rotate a b c]
  simp only [hα, hβ, hγ]
  field_simp
  ring

namespace TriMesh

open Matrix

variable {V F : Type*} (M : TriMesh V F)

/-- The discrete conformal change of metric with log scale factors `u`:
`ℓ̃_ij = e^{(u_i + u_j)/2}·ℓ_ij` (the `u` of `CLAUDE.md`). -/
def conformalScale (u : V → ℝ) : TriMesh V F :=
  { M with len := fun i j => M.len i j * exp ((u i + u j) / 2) }

/-- Conformal changes compose additively in `u`. -/
theorem conformalScale_conformalScale (u w : V → ℝ) :
    (M.conformalScale u).conformalScale w = M.conformalScale (fun i => u i + w i) := by
  simp only [conformalScale, mk.injEq, true_and]
  funext i j
  rw [mul_assoc, ← exp_add]
  congr 2
  ring

/-- The zero scale factor does not change the mesh. -/
theorem conformalScale_zero : M.conformalScale (fun _ => 0) = M := by
  simp [conformalScale]

variable [DecidableEq V]

/-- Row `i` of the weights of one face applied to a vector. -/
theorem sum_cotOpp_mul [Fintype V] (hM : M.Valid) (f : F) (i : V) (g : V → ℝ) :
    ∑ j, M.cotOpp f i j * g j =
      M.cotC f * ((if i = M.a f then 1 else 0) * g (M.b f) +
        (if i = M.b f then 1 else 0) * g (M.a f)) +
      M.cotA f * ((if i = M.b f then 1 else 0) * g (M.c f) +
        (if i = M.c f then 1 else 0) * g (M.b f)) +
      M.cotB f * ((if i = M.c f then 1 else 0) * g (M.a f) +
        (if i = M.a f then 1 else 0) * g (M.c f)) := by
  have key : ∀ j, M.cotOpp f i j * g j =
      M.cotC f * ((if i = M.a f then (1 : ℝ) else 0) * ((if j = M.b f then (1 : ℝ) else 0) * g j)) +
      M.cotC f * ((if i = M.b f then (1 : ℝ) else 0) * ((if j = M.a f then (1 : ℝ) else 0) * g j)) +
      M.cotA f * ((if i = M.b f then (1 : ℝ) else 0) * ((if j = M.c f then (1 : ℝ) else 0) * g j)) +
      M.cotA f * ((if i = M.c f then (1 : ℝ) else 0) * ((if j = M.b f then (1 : ℝ) else 0) * g j)) +
      M.cotB f * ((if i = M.c f then (1 : ℝ) else 0) * ((if j = M.a f then (1 : ℝ) else 0) * g j)) +
      M.cotB f * ((if i = M.a f then (1 : ℝ) else 0) * ((if j = M.c f then (1 : ℝ) else 0) * g j)) := by
    intro j
    rw [M.cotOpp_eq hM]
    ring
  have hind : ∀ y : V, ∑ j, (if j = y then (1 : ℝ) else 0) * g j = g y := fun y => by simp
  simp only [key, Finset.sum_add_distrib, ← Finset.mul_sum, hind]
  ring

/-- First variation of the corner angle of the face `f` at the vertex `i` under the conformal
change `t·v`: `½ Σ_j cotOpp f i j·(v_j − v_i)`. -/
theorem hasDerivAt_angleAt_conformalScale [Fintype V] (hM : M.Valid) (v : V → ℝ) (f : F) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleAt f i)
      (1 / 2 * ∑ j, M.cotOpp f i j * (v j - v i)) 0 := by
  have hab := hM.ab f
  have hbc := hM.bc f
  have hca := hM.ca f
  have ht := hM.tri f
  rw [M.sum_cotOpp_mul hM]
  have eA : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleA f) =
      fun t : ℝ => cornerAngle (M.la f * exp (t * ((v (M.b f) + v (M.c f)) / 2)))
        (M.lb f * exp (t * ((v (M.c f) + v (M.a f)) / 2)))
        (M.lc f * exp (t * ((v (M.a f) + v (M.b f)) / 2))) := by
    funext t
    simp only [angleA, la, lb, lc, conformalScale]
    congr 3 <;> ring
  have eB : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleB f) =
      fun t : ℝ => cornerAngle (M.lb f * exp (t * ((v (M.c f) + v (M.a f)) / 2)))
        (M.lc f * exp (t * ((v (M.a f) + v (M.b f)) / 2)))
        (M.la f * exp (t * ((v (M.b f) + v (M.c f)) / 2))) := by
    funext t
    simp only [angleB, la, lb, lc, conformalScale]
    congr 3 <;> ring
  have eC : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleC f) =
      fun t : ℝ => cornerAngle (M.lc f * exp (t * ((v (M.a f) + v (M.b f)) / 2)))
        (M.la f * exp (t * ((v (M.b f) + v (M.c f)) / 2)))
        (M.lb f * exp (t * ((v (M.c f) + v (M.a f)) / 2))) := by
    funext t
    simp only [angleC, la, lb, lc, conformalScale]
    congr 3 <;> ring
  by_cases hia : i = M.a f
  · have hfun : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleAt f i) =
        fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleA f := by
      funext t
      simp [angleAt, conformalScale, hia]
    rw [hfun, eA]
    refine (hasDerivAt_cornerAngle_conformal ht _ _ _).congr_deriv ?_
    simp only [hia, cotA, cotB, cotC, hab, ite_true, ite_false, hca.symm]
    ring
  · by_cases hib : i = M.b f
    · have hfun : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleAt f i) =
          fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleB f := by
        funext t
        simp [angleAt, conformalScale, hib, hab.symm]
      rw [hfun, eB]
      refine (hasDerivAt_cornerAngle_conformal ht.rotate _ _ _).congr_deriv ?_
      simp only [hib, cotA, cotB, cotC, hab.symm, ite_true, ite_false, hbc]
      ring
    · by_cases hic : i = M.c f
      · have hfun : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleAt f i) =
            fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleC f := by
          funext t
          simp [angleAt, conformalScale, hic, hca, hbc.symm]
        rw [hfun, eC]
        refine (hasDerivAt_cornerAngle_conformal ht.rotate.rotate _ _ _).congr_deriv ?_
        simp only [hic, cotA, cotB, cotC, hca, hbc.symm, ite_true, ite_false]
        ring
      · have hfun : (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleAt f i) =
            fun _ => (0 : ℝ) := by
          funext t
          simp [angleAt, conformalScale, hia, hib, hic]
        rw [hfun]
        simp only [hia, hib, hic, ite_false, zero_mul, add_zero, mul_zero]
        exact hasDerivAt_const _ _

variable [Fintype V] [Fintype F]

/-- **The linearized discrete Yamabe equation.** Under the conformal change
`ℓ̃_ij = e^{t(v_i + v_j)/2}·ℓ_ij`, the angle sum at the vertex `i` has derivative `−(L v)_i` at
`t = 0`, where `L` is the cotan Laplacian: `dΘ = −L du`.

[formula: discrete-yamabe-linearization] -/
theorem hasDerivAt_angleSum_conformalScale (hM : M.Valid) (v : V → ℝ) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleSum i)
      (-(M.laplacian *ᵥ v) i) 0 := by
  have hsum : HasDerivAt
      (fun t : ℝ => ∑ f, (M.conformalScale (fun j => t * v j)).angleAt f i)
      (∑ f, 1 / 2 * ∑ j, M.cotOpp f i j * (v j - v i)) 0 :=
    HasDerivAt.fun_sum fun f _ => M.hasDerivAt_angleAt_conformalScale hM v f i
  refine hsum.congr_deriv ?_
  rw [laplacian_mulVec_apply, ← Finset.sum_neg_distrib]
  simp only [offDiag, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [neg_div, neg_neg, Finset.sum_div, Finset.sum_mul, ← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl fun f _ => by ring

/-- **Linearized discrete Yamabe equation, curvature form.** The angle defect of the
conformally scaled metric has derivative `(L v)_i`: `K̃(t v) = Ω + t·L v + o(t)`. Asking the
scaled metric to have curvature `μ` therefore gives, to first order, the state equation
`L u = μ − Ω` of the canonical problem.

[formula: discrete-yamabe-linearization] -/
theorem hasDerivAt_angleDefect_conformalScale (hM : M.Valid) (v : V → ℝ) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => t * v j)).angleDefect i)
      ((M.laplacian *ᵥ v) i) 0 := by
  have h := (M.hasDerivAt_angleSum_conformalScale hM v i).const_sub (2 * π)
  simpa [angleDefect] using h

/-- The same at an arbitrary scale factor `u`: the derivative of the angle defect of
`e^{(u_i + u_j)/2} ℓ_ij` in the direction `v` is `(L̃ v)_i`, with `L̃` the cotan Laplacian of the
scaled metric (assumed to satisfy the triangle inequalities).

[formula: discrete-yamabe-linearization] -/
theorem hasDerivAt_angleDefect_conformalScale_add (u : V → ℝ) (hM : (M.conformalScale u).Valid)
    (v : V → ℝ) (i : V) :
    HasDerivAt (fun t : ℝ => (M.conformalScale (fun j => u j + t * v j)).angleDefect i)
      (((M.conformalScale u).laplacian *ᵥ v) i) 0 := by
  have h := (M.conformalScale u).hasDerivAt_angleDefect_conformalScale hM v i
  simpa only [conformalScale_conformalScale] using h

end TriMesh

end DDG
