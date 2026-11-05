import DdgProofs.RepulsiveCurves.Hierarchy

/-!
# Repulsive Curves: the discrete adjoint of the Barnes–Hut energy (SPEC §17i)

On fixed interaction lists the evaluated Barnes–Hut energy is a sum of exact pair terms `e_IJ`
(whose derivative is `rc-energy-differential`, `Energy.lean`) and lumped terms
`ℓ_I L_N k(x_I, x̄_N, T_I)` (Eq. 26, `lumpedTerm` of `Hierarchy.lean`) whose node data
`L_N = Σ_{J∈N} ℓ_J`, `x̄_N = Σ ℓ_J x_J / L_N` move with the edges of the node.

* `hasDerivAt_nodeCentre`: `x̄_N' = L_N⁻¹ Σ_{J∈N} (ℓ_J' (x_J − x̄_N) + ℓ_J x_J')`.
* `bh_adjoint` (`rc-bh-adjoint`): the derivative of one lumped term along a motion of the curve
  is `Σ_a g_a · Γ'_{i_a} + Σ_{J∈N} Σ_b h_{J,b} · Γ'_{j_b}` with the own part
  `g_a = L_N (s_a (k T_I + ∇_T k) + (ℓ_I/2) ∇_p k)` and, for the edges of the node,
  `h_{J,b} = s_b λ_J T_J + (ℓ_J/2) u`, `u = −ℓ_I ∇_p k`, `λ_J = ℓ_I k + u · (x_J − x̄_N)`
  (all kernel quantities at `(x_I, x̄_N, T_I)`).
* `bh_adjoint_shift` (`rc-bh-adjoint-shift`): the downward pass. Along the path from the root to
  a leaf, `β ← β + U · (c' − c) + b`, `U ← U + u` (re-centred at each node's centre `c`) gives
  `β + U · (x − c_leaf) = Σ_nodes (b + u · (x − c))` and `U = Σ u`; so the per-node sums of the
  lumped terms' `(b, u)` reach every edge of the node with its own `λ_J` and `μ_J = U`.

Summing over the lumped terms of all rows and over the exact pairs (sum rule) gives the
derivative of the evaluated energy on the lists; that assembly is linear and not restated here.
-/

noncomputable section

namespace DDG

namespace RC

open Matrix Real Finset

/-! ## Derivatives of the node data -/

section nodeDeriv

variable {E : Type*} {N : Finset E} {ℓ : ℝ → E → ℝ} {x : ℝ → E → V3} {ℓ' : E → ℝ}
  {x' : E → V3} {t : ℝ}

/-- `L_N' = Σ_{J∈N} ℓ_J'`.

[formula: rc-bh-adjoint] -/
theorem hasDerivAt_nodeMass (hℓ : ∀ J ∈ N, HasDerivAt (fun s => ℓ s J) (ℓ' J) t) :
    HasDerivAt (fun s => nodeMass (ℓ s) N) (∑ J ∈ N, ℓ' J) t := by
  unfold nodeMass
  exact HasDerivAt.fun_sum (A := fun J s => ℓ s J) hℓ

/-- `(Σ ℓ_J x_J)' = Σ (ℓ_J x_J' + ℓ_J' x_J)`.

[formula: rc-bh-adjoint] -/
theorem hasDerivAt_nodeMoment (hℓ : ∀ J ∈ N, HasDerivAt (fun s => ℓ s J) (ℓ' J) t)
    (hx : ∀ J ∈ N, HasDerivAt (fun s => x s J) (x' J) t) :
    HasDerivAt (fun s => nodeMoment (ℓ s) (x s) N)
      (∑ J ∈ N, (ℓ t J • x' J + ℓ' J • x t J)) t := by
  unfold nodeMoment
  exact HasDerivAt.fun_sum (A := fun J s => ℓ s J • x s J) fun J hJ => (hℓ J hJ).smul (hx J hJ)

/-- **Derivative of the node centre**: `x̄_N' = L_N⁻¹ Σ_{J∈N} (ℓ_J' (x_J − x̄_N) + ℓ_J x_J')`
(`L_N ≠ 0`).

[formula: rc-bh-adjoint] -/
theorem hasDerivAt_nodeCentre (hℓ : ∀ J ∈ N, HasDerivAt (fun s => ℓ s J) (ℓ' J) t)
    (hx : ∀ J ∈ N, HasDerivAt (fun s => x s J) (x' J) t) (hL : nodeMass (ℓ t) N ≠ 0) :
    HasDerivAt (fun s => nodeCentre (ℓ s) (x s) N)
      ((nodeMass (ℓ t) N)⁻¹ •
        ∑ J ∈ N, (ℓ' J • (x t J - nodeCentre (ℓ t) (x t) N) + ℓ t J • x' J)) t := by
  have h := ((hasDerivAt_nodeMass hℓ).inv hL).smul (hasDerivAt_nodeMoment hℓ hx)
  simp only [Pi.inv_apply] at h
  unfold nodeCentre
  refine h.congr_deriv ?_
  set L := nodeMass (ℓ t) N
  set M := nodeMoment (ℓ t) (x t) N
  have hsum : ∑ J ∈ N, (ℓ' J • (x t J - L⁻¹ • M) + ℓ t J • x' J) =
      ∑ J ∈ N, (ℓ t J • x' J + ℓ' J • x t J) - (∑ J ∈ N, ℓ' J) • (L⁻¹ • M) := by
    rw [Finset.sum_smul, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun J _ => ?_
    rw [smul_sub]
    abel
  rw [hsum, smul_sub, smul_smul, smul_smul]
  ext i
  simp only [Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  field_simp
  ring

end nodeDeriv

/-! ## Two algebraic identities -/

section alg

/-- `∇_T k ⊥ T`. -/
theorem kgradT_dot_tangent {α β : ℝ} {p q T : V3} (hT : T ⬝ᵥ T = 1) :
    kgradT α β p q T ⬝ᵥ T = 0 := by
  rw [kgradT, smul_dotProduct, dotProduct_comm (proj T (p - q)) T, dot_proj_self hT, smul_zero]

/-- The own part: with `x_I' = ½(g₀ + g₁)`, `T_I' = ℓ_I⁻¹ Σ_a s_a P_T g_a`, `ℓ_I' = Σ_a s_a T·g_a`
and `∇_T k ⊥ T`, `(∇_p k · x_I' + ∇_T k · T_I') ℓ_I L + k ℓ_I' L = Σ_a g_a · (own part)_a`.

[formula: rc-bh-adjoint] -/
theorem adj_own_alg {lI L k : ℝ} {TI gP gT : V3} (g : Fin 2 → V3) (hl0 : lI ≠ 0)
    (hgT : gT ⬝ᵥ TI = 0) :
    (gP ⬝ᵥ ((1 / 2 : ℝ) • (g 0 + g 1)) + gT ⬝ᵥ (lI⁻¹ • ∑ a, sgn a • proj TI (g a))) * lI * L
        + k * (∑ a, sgn a * (TI ⬝ᵥ g a)) * L =
      ∑ a, (L • (sgn a • (k • TI + gT) + (lI / 2) • gP)) ⬝ᵥ g a := by
  simp only [Fin.sum_univ_two, sgn_zero, sgn_one, dotProduct_smul, dotProduct_add, smul_dotProduct,
    add_dotProduct, dot_proj, hgT, smul_eq_mul]
  field_simp
  ring

/-- The node part: with `x̄' = L⁻¹ Σ_J (ℓ_J' (x_J − c) + ℓ_J x_J')`, `ℓ_J' = Σ_b s_b T_J·g_{J,b}`,
`x_J' = ½(g_{J,0} + g_{J,1})` and `L' = Σ_J ℓ_J'`,
`−(∇_p k · x̄') ℓ_I L + k ℓ_I L' = Σ_J Σ_b (s_b λ_J T_J + (ℓ_J/2) u) · g_{J,b}` with
`u = −ℓ_I ∇_p k`, `λ_J = ℓ_I k + u·(x_J − c)`.

[formula: rc-bh-adjoint] -/
theorem adj_node_alg {E : Type*} {L lI k : ℝ} {gP c : V3} (hL : L ≠ 0) (N : Finset E)
    (T x : E → V3) (l : E → ℝ) (g : E → Fin 2 → V3) :
    -(gP ⬝ᵥ (L⁻¹ • ∑ J ∈ N, ((∑ a, sgn a * (T J ⬝ᵥ g J a)) • (x J - c)
          + l J • ((1 / 2 : ℝ) • (g J 0 + g J 1))))) * lI * L
        + k * lI * ∑ J ∈ N, ∑ a, sgn a * (T J ⬝ᵥ g J a) =
      ∑ J ∈ N, ∑ b, ((sgn b * (lI * k + ((-lI) • gP) ⬝ᵥ (x J - c))) • T J
        + (l J / 2) • ((-lI) • gP)) ⬝ᵥ g J b := by
  rw [dotProduct_smul, dotProduct_sum, smul_eq_mul, Finset.mul_sum, Finset.mul_sum,
    ← Finset.sum_neg_distrib, Finset.sum_mul, Finset.sum_mul, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun J _ => ?_
  simp only [Fin.sum_univ_two, sgn_zero, sgn_one, dotProduct_smul, dotProduct_add, smul_dotProduct,
    add_dotProduct, smul_eq_mul, neg_mul, one_mul]
  field_simp
  ring

end alg

/-! ## The adjoint of one lumped term -/

section adjoint

variable {V E : Type*} {Γ : ℝ → V → V3} {Γ' : V → V3} {t : ℝ} {ep : E → Fin 2 → V}

/-- The own part `g_a = L_N (s_a (k T_I + ∇_T k) + (ℓ_I/2) ∇_p k)` at `(x_I, x̄, T_I)`.

[formula: rc-bh-adjoint] -/
def adjOwn (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I : E) (L : ℝ) (c : V3) (a : Fin 2) :
    V3 :=
  L • (sgn a • (kernel α β (mid γ ep I) c (tangent γ ep I) • tangent γ ep I
      + kgradT α β (mid γ ep I) c (tangent γ ep I))
    + (len γ ep I / 2) • kgradP α β (mid γ ep I) c (tangent γ ep I))

/-- `u = −ℓ_I ∇_p k(x_I, x̄, T_I)`: the derivative of the lumped term with respect to `x̄`
divided by `L_N`.

[formula: rc-bh-adjoint] -/
def adjU (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I : E) (c : V3) : V3 :=
  (-len γ ep I) • kgradP α β (mid γ ep I) c (tangent γ ep I)

/-- `λ_J = ℓ_I k + u · (x_J − x̄)`: the derivative with respect to `ℓ_J`.

[formula: rc-bh-adjoint] -/
def adjLambda (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) (c : V3) : ℝ :=
  len γ ep I * kernel α β (mid γ ep I) c (tangent γ ep I) + adjU α β γ ep I c ⬝ᵥ (mid γ ep J - c)

/-- The node part `h_{J,b} = s_b λ_J T_J + (ℓ_J/2) u`.

[formula: rc-bh-adjoint] -/
def adjNode (α β : ℝ) (γ : V → V3) (ep : E → Fin 2 → V) (I J : E) (c : V3) (b : Fin 2) : V3 :=
  (sgn b * adjLambda α β γ ep I J c) • tangent γ ep J + (len γ ep J / 2) • adjU α β γ ep I c

/-- **The derivative of a lumped Barnes–Hut term with moving node data** (SPEC §17i): along a
differentiable motion `Γ` with `ℓ_I > 0`, `ℓ_J > 0` for the edges of `N`, `L_N ≠ 0` and
`ν = |T_I × (x_I − x̄_N)| > 0`,
`d/ds [ℓ_I L_N k(x_I, x̄_N, T_I)] = Σ_a g_a · Γ'_{i_a} + Σ_{J∈N} Σ_b h_{J,b} · Γ'_{j_b}`.

[formula: rc-bh-adjoint] -/
theorem bh_adjoint {α β : ℝ} (hΓ : ∀ v, HasDerivAt (fun s => Γ s v) (Γ' v) t) {I : E}
    {N : Finset E} (hI : 0 < len (Γ t) ep I) (hN : ∀ J ∈ N, 0 < len (Γ t) ep J)
    (hL : nodeMass (len (Γ t) ep) N ≠ 0)
    (hν : crossProduct (tangent (Γ t) ep I)
      (mid (Γ t) ep I - nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N) ≠ 0) :
    HasDerivAt (fun s => lumpedTerm (len (Γ s) ep) (mid (Γ s) ep) α β (mid (Γ s) ep I)
        (tangent (Γ s) ep I) (len (Γ s) ep I) N)
      (∑ a, adjOwn α β (Γ t) ep I (nodeMass (len (Γ t) ep) N)
          (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N) a ⬝ᵥ Γ' (ep I a)
        + ∑ J ∈ N, ∑ b, adjNode α β (Γ t) ep I J
          (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N) b ⬝ᵥ Γ' (ep J b)) t := by
  -- Edge derivatives.
  have hlen : ∀ J ∈ N, HasDerivAt (fun s => len (Γ s) ep J)
      (∑ a, sgn a * (tangent (Γ t) ep J ⬝ᵥ Γ' (ep J a))) t :=
    fun J hJ => hasDerivAt_len (ep := ep) (I := J) hΓ (hN J hJ)
  have hmid : ∀ J ∈ N, HasDerivAt (fun s => mid (Γ s) ep J)
      ((1 / 2 : ℝ) • (Γ' (ep J 0) + Γ' (ep J 1))) t :=
    fun J _ => hasDerivAt_mid (ep := ep) (I := J) hΓ
  have hc := hasDerivAt_nodeCentre (ℓ := fun s => len (Γ s) ep) (x := fun s => mid (Γ s) ep)
    hlen hmid hL
  have hM := hasDerivAt_nodeMass (ℓ := fun s => len (Γ s) ep) hlen
  have hT := hasDerivAt_tangent (ep := ep) (I := I) hΓ hI
  have hx := hasDerivAt_mid (ep := ep) (I := I) hΓ
  have hunit := tangent_dot_self hI
  have hk := hasDerivAt_kernel (α := α) (β := β) hx hc hT hunit (tangent_dot_deriv hI) hν
  have hl := hasDerivAt_len (ep := ep) (I := I) hΓ hI
  unfold lumpedTerm
  refine ((hk.mul hl).mul hM).congr_deriv ?_
  have eqI := adj_own_alg (lI := len (Γ t) ep I) (L := nodeMass (len (Γ t) ep) N)
    (k := kernel α β (mid (Γ t) ep I) (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N)
      (tangent (Γ t) ep I))
    (gP := kgradP α β (mid (Γ t) ep I) (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N)
      (tangent (Γ t) ep I))
    (gT := kgradT α β (mid (Γ t) ep I) (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N)
      (tangent (Γ t) ep I))
    (fun a => Γ' (ep I a)) hI.ne' (kgradT_dot_tangent hunit)
  have eqN := adj_node_alg (lI := len (Γ t) ep I)
    (k := kernel α β (mid (Γ t) ep I) (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N)
      (tangent (Γ t) ep I))
    (gP := kgradP α β (mid (Γ t) ep I) (nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N)
      (tangent (Γ t) ep I))
    (c := nodeCentre (len (Γ t) ep) (mid (Γ t) ep) N) hL N (fun J => tangent (Γ t) ep J)
    (fun J => mid (Γ t) ep J) (fun J => len (Γ t) ep J) (fun J b => Γ' (ep J b))
  beta_reduce at eqI eqN
  simp only [adjOwn, adjNode, adjLambda, adjU, Pi.mul_apply]
  rw [dotProduct_sub]
  linear_combination eqI + eqN

end adjoint

/-! ## The downward pass -/

section shift

/-- One step of the downward pass: the state `(β, U, c)` (accumulated scalar, accumulated vector,
current centre) moves to a child node with data `(b, u, c')`:
`(β + U · (c' − c) + b, U + u, c')`.

[formula: rc-bh-adjoint-shift] -/
def shiftStep (st : ℝ × V3 × V3) (nd : ℝ × V3 × V3) : ℝ × V3 × V3 :=
  (st.1 + st.2.1 ⬝ᵥ (nd.2.2 - st.2.2) + nd.1, st.2.1 + nd.2.1, nd.2.2)

/-- **The re-centred downward pass sums the node terms** (SPEC §17i): folding `shiftStep` along
a path of nodes `(b_k, u_k, c_k)` from the state `(β₀, U₀, c₀)` gives `(β, U, c)` with
`β + U · (x − c) = β₀ + U₀ · (x − c₀) + Σ_k (b_k + u_k · (x − c_k))` and `U = U₀ + Σ_k u_k`
for every point `x`. From the state `(0, 0, c₀)` at the root this is the sum over the path.

[formula: rc-bh-adjoint-shift] -/
theorem bh_adjoint_shift (path : List (ℝ × V3 × V3)) (st : ℝ × V3 × V3) (x : V3) :
    (path.foldl shiftStep st).1 + (path.foldl shiftStep st).2.1 ⬝ᵥ
        (x - (path.foldl shiftStep st).2.2) =
      st.1 + st.2.1 ⬝ᵥ (x - st.2.2) + (path.map fun nd => nd.1 + nd.2.1 ⬝ᵥ (x - nd.2.2)).sum ∧
    (path.foldl shiftStep st).2.1 = st.2.1 + (path.map fun nd => nd.2.1).sum := by
  induction path generalizing st with
  | nil => simp
  | cons nd rest ih =>
    obtain ⟨h1, h2⟩ := ih (shiftStep st nd)
    simp only [List.foldl_cons, List.map_cons, List.sum_cons]
    refine ⟨?_, ?_⟩
    · rw [h1]
      simp only [shiftStep, add_dotProduct, dotProduct_sub]
      ring
    · rw [h2]
      simp only [shiftStep]
      abel

/-- From the zero state at the root: `β + U · (x − c_leaf) = Σ_path (b + u · (x − c))`.

[formula: rc-bh-adjoint-shift] -/
theorem bh_adjoint_shift_root (path : List (ℝ × V3 × V3)) (c₀ x : V3) :
    (path.foldl shiftStep (0, 0, c₀)).1 + (path.foldl shiftStep (0, 0, c₀)).2.1 ⬝ᵥ
        (x - (path.foldl shiftStep (0, 0, c₀)).2.2) =
      (path.map fun nd => nd.1 + nd.2.1 ⬝ᵥ (x - nd.2.2)).sum ∧
    (path.foldl shiftStep (0, 0, c₀)).2.1 = (path.map fun nd => nd.2.1).sum := by
  obtain ⟨h1, h2⟩ := bh_adjoint_shift path (0, 0, c₀) x
  simp only [zero_dotProduct, zero_add] at h1 h2
  exact ⟨h1, h2⟩

end shift

end RC

end DDG
