import DDG.FEM.Elements
import DDG.FEM.Transfer

/-!
# Integrals over a tiling by triangles

A finite family of nondegenerate triangles `T_k = triangle (q k)`, `k ∈ s`, **tiles** a set `Ω` of
the plane when their union is `Ω` and their interiors are pairwise disjoint (`IsTiling`). That is
the hypothesis; nothing here proves that a given mesh layout tiles a domain.

* `volume_frontier_triangle` — the boundary of a triangle is a null set (it is the frontier of a
  convex set, `Convex.addHaar_frontier`); hence two triangles with disjoint interiors meet in a
  null set (`aedisjoint_triangle`);
* `integral_biUnion_triangle`, `IsTiling.integral_eq_sum` — **additivity**:
  `∫_Ω f = Σ_k ∫_{T_k} f` for `f` integrable on every triangle (in particular for continuous `f`,
  `IsTiling.integral_eq_sum_of_continuous`), and `|Ω| = Σ_k |T_k|` (`IsTiling.volume_real_eq_sum`).

On it:

* **`ic-transfer-mass-integral`** — a piece `r` of the common subdivision inside the input triangle
  `p` has the barycentric corner rows `β_{j c} = λ_c(r_j)`; its signed double area is
  `det β · D(p)` (`D_eq_bdet_mul`), the hat functions of `p` restricted to it are the linear
  interpolants of the rows (`bary_eq_interp`), and so
  `∫_r λ_a λ_b = |det β| Area(p)/12 · Σ_{j,l} (1 + δ_jl) β_{ja} β_{lb}` (`integral_piece_mass`),
  the contribution of the piece to `(P₀ᵀ M_S P₀)_{ab}` as `Transfer::new` assembles it (fan
  triangle area `|det β| · Area(T)`, entries `A/12 (1 + δ)`, rows of `P₀` the `β`). Summed over
  pieces that tile `p` this is `∫_p λ_a λ_b = Area(p)(1 + δ_ab)/12`, the entry of `M₀`
  (`transfer_mass_tiling`); for pieces of the orientation of `p`, `Σ pieceMass β a b = (1 + δ_ab)/12`
  (`sum_pieceMass_tiling`).

Chen–Holst's interpolation error over a tiling (`it-odt-energy`) is in `DDG/Delaunay/ODT.lean`.
-/

noncomputable section

namespace DDG

namespace FEM

open MeasureTheory Set Real

/-! ## The boundary of a triangle is null -/

theorem convex_triangle (p : Fin 3 → ℝ × ℝ) : Convex ℝ (triangle p) := by
  rw [triangle_eq_convexHull]
  exact convex_convexHull ℝ _

theorem isClosed_triangle (p : Fin 3 → ℝ × ℝ) : IsClosed (triangle p) :=
  (isCompact_triangle p).isClosed

/-- **The boundary of a triangle has measure zero.** -/
theorem volume_frontier_triangle (p : Fin 3 → ℝ × ℝ) : volume (frontier (triangle p)) = 0 :=
  (convex_triangle p).addHaar_frontier volume

/-- Two triangles with disjoint interiors meet on their boundaries. -/
theorem triangle_inter_subset {p q : Fin 3 → ℝ × ℝ}
    (h : Disjoint (interior (triangle p)) (interior (triangle q))) :
    triangle p ∩ triangle q ⊆ frontier (triangle p) ∪ frontier (triangle q) := by
  rintro x ⟨hp, hq⟩
  by_cases hxp : x ∈ interior (triangle p)
  · exact Or.inr ⟨subset_closure hq, fun hxq => Set.disjoint_left.1 h hxp hxq⟩
  · exact Or.inl ⟨subset_closure hp, hxp⟩

/-- **Two triangles with disjoint interiors intersect in a null set.** -/
theorem aedisjoint_triangle {p q : Fin 3 → ℝ × ℝ}
    (h : Disjoint (interior (triangle p)) (interior (triangle q))) :
    AEDisjoint volume (triangle p) (triangle q) :=
  measure_mono_null (triangle_inter_subset h)
    (measure_union_null (volume_frontier_triangle p) (volume_frontier_triangle q))

/-! ## Additivity -/

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- **Additivity of the integral over triangles with pairwise disjoint interiors.** -/
theorem integral_biUnion_triangle {ι : Type*} (s : Finset ι) (q : ι → Fin 3 → ℝ × ℝ)
    (hd : (s : Set ι).Pairwise fun k l => Disjoint (interior (triangle (q k))) (interior (triangle (q l))))
    {f : ℝ × ℝ → E} (hf : ∀ k ∈ s, IntegrableOn f (triangle (q k))) :
    ∫ x in ⋃ k ∈ s, triangle (q k), f x = ∑ k ∈ s, ∫ x in triangle (q k), f x := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha IH =>
    simp only [Finset.coe_insert, Set.pairwise_insert, Finset.set_biUnion_insert,
      Finset.forall_mem_insert] at hd hf ⊢
    rw [setIntegral_union₀ _ _ hf.1 (integrableOn_finset_iUnion.2 hf.2), Finset.sum_insert ha,
      IH hd.1 hf.2]
    · show volume (triangle (q a) ∩ ⋃ i ∈ s, triangle (q i)) = 0
      rw [Set.inter_iUnion₂]
      refine le_antisymm ((measure_biUnion_finset_le s _).trans ?_) (zero_le)
      refine (Finset.sum_eq_zero fun i hi => ?_).le
      exact aedisjoint_triangle (hd.2 i (Finset.mem_coe.2 hi) (fun h => ha (h ▸ hi))).1
    · exact (Finset.measurableSet_biUnion _ fun i _ => measurableSet_triangle _).nullMeasurableSet

/-- A finite family of nondegenerate triangles `T_k = triangle (q k)` (`k ∈ s`) **tiles** `Ω`:
their union is `Ω` and their interiors are pairwise disjoint (either orientation). -/
structure IsTiling {ι : Type*} (s : Finset ι) (q : ι → Fin 3 → ℝ × ℝ) (Ω : Set (ℝ × ℝ)) : Prop where
  nondeg : ∀ k ∈ s, D (q k) ≠ 0
  cover : (⋃ k ∈ s, triangle (q k)) = Ω
  disjoint : (s : Set ι).Pairwise fun k l =>
    Disjoint (interior (triangle (q k))) (interior (triangle (q l)))

namespace IsTiling

variable {ι : Type*} {s : Finset ι} {q : ι → Fin 3 → ℝ × ℝ} {Ω : Set (ℝ × ℝ)}

theorem triangle_subset (h : IsTiling s q Ω) {k : ι} (hk : k ∈ s) : triangle (q k) ⊆ Ω := by
  rw [← h.cover]
  exact Set.subset_biUnion_of_mem (u := fun k => triangle (q k)) hk

theorem measurableSet (h : IsTiling s q Ω) : MeasurableSet Ω := by
  rw [← h.cover]
  exact Finset.measurableSet_biUnion _ fun i _ => measurableSet_triangle _

/-- **The integral over a tiled set is the sum of the integrals over the tiles**, for `f`
integrable on every tile. -/
theorem integral_eq_sum (h : IsTiling s q Ω) {f : ℝ × ℝ → E}
    (hf : ∀ k ∈ s, IntegrableOn f (triangle (q k))) :
    ∫ x in Ω, f x = ∑ k ∈ s, ∫ x in triangle (q k), f x := by
  rw [← h.cover]
  exact integral_biUnion_triangle s q h.disjoint hf

theorem integral_eq_sum_of_continuous (h : IsTiling s q Ω) {f : ℝ × ℝ → E} (hf : Continuous f) :
    ∫ x in Ω, f x = ∑ k ∈ s, ∫ x in triangle (q k), f x :=
  h.integral_eq_sum fun _ _ => integrableOn_triangle _ hf

omit [NormedSpace ℝ E] in
theorem integrableOn (h : IsTiling s q Ω) {f : ℝ × ℝ → E} (hf : Continuous f) :
    IntegrableOn f Ω := by
  rw [← h.cover]
  exact integrableOn_finset_iUnion.2 fun _ _ => integrableOn_triangle _ hf

/-- **The area of a tiled set is the sum of the areas of the tiles.** -/
theorem volume_real_eq_sum (h : IsTiling s q Ω) : volume.real Ω = ∑ k ∈ s, area (q k) := by
  have e := h.integral_eq_sum_of_continuous (f := fun _ => (1 : ℝ)) continuous_const
  rw [setIntegral_const, smul_eq_mul, mul_one] at e
  rw [e]
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [setIntegral_const, smul_eq_mul, mul_one, volume_real_triangle_eq_area (h.nondeg k hk)]

end IsTiling

/-! ## `ic-transfer-mass-integral`: the transfer mass matrix as an integral -/

section Transfer


/-- The signed double area of three affine combinations of the corners: `det β · D(p)`. -/
theorem sda_affine (p : Fin 3 → ℝ × ℝ) (β : Fin 3 → Fin 3 → ℝ)
    (hβ : ∀ j, β j 0 + β j 1 + β j 2 = 1) :
    sda (∑ c, β 0 c • p c) (∑ c, β 1 c • p c) (∑ c, β 2 c • p c) = bdet β * D p := by
  have e0 : β 0 2 = 1 - β 0 0 - β 0 1 := by linarith [hβ 0]
  have e1 : β 1 2 = 1 - β 1 0 - β 1 1 := by linarith [hβ 1]
  have e2 : β 2 2 = 1 - β 2 0 - β 2 1 := by linarith [hβ 2]
  simp only [sda, D, bdet, Fin.sum_univ_three, Prod.fst_add, Prod.snd_add, Prod.smul_fst,
    Prod.smul_snd, smul_eq_mul, e0, e1, e2]
  ring

/-- The barycentric rows of a piece `r` in the triangle `p`: `β_{j c} = λ_c(r_j)`. -/
def pieceRows (p r : Fin 3 → ℝ × ℝ) : Fin 3 → Fin 3 → ℝ := fun j c => bary p c (r j)

/-- **The signed double area of a piece is `det β · D(p)`** (so its area is `|det β| · Area(p)`, as
`triangle_areas` in `transfer.rs` computes it). -/
theorem D_eq_bdet_mul {p : Fin 3 → ℝ × ℝ} (hp : D p ≠ 0) (r : Fin 3 → ℝ × ℝ) :
    D r = bdet (pieceRows p r) * D p := by
  have hr : ∀ j, r j = ∑ c, pieceRows p r j c • p c := fun j => by
    simp only [pieceRows, Fin.sum_univ_three]
    exact (sum_bary_smul (r j) hp).symm
  have := sda_affine p (pieceRows p r) fun j => sum_bary (r j) hp
  rw [← hr 0, ← hr 1, ← hr 2] at this
  exact this

theorem area_piece {p : Fin 3 → ℝ × ℝ} (hp : D p ≠ 0) (r : Fin 3 → ℝ × ℝ) :
    area r = |bdet (pieceRows p r)| * area p := by
  have h1 := abs_D_eq_two_mul_area r
  have h2 := abs_D_eq_two_mul_area p
  rw [D_eq_bdet_mul hp r, abs_mul] at h1
  linarith [h1, h2, show |bdet (pieceRows p r)| * |D p| = |bdet (pieceRows p r)| * (2 * area p)
    by rw [h2]]

/-- **On a piece, the hat functions of `p` are the linear interpolants of their corner values**
(`P₀` represents them exactly). -/
theorem bary_eq_interp (p : Fin 3 → ℝ × ℝ) {r : Fin 3 → ℝ × ℝ} (hr : D r ≠ 0) (a : Fin 3)
    (x : ℝ × ℝ) : bary p a x = interp r (fun j => bary p a (r j)) x := by
  have hs := sum_bary x hr
  have hx := sum_bary_smul x hr
  have hx1 := congrArg Prod.fst hx
  have hx2 := congrArg Prod.snd hx
  simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] at hx1 hx2
  simp only [interp, Fin.sum_univ_three]
  rw [bary_affine p a x, bary_affine p a (r 0), bary_affine p a (r 1), bary_affine p a (r 2)]
  linear_combination (-(bary p a 0)) * hs - (baryGrad p a).1 * hx1 - (baryGrad p a).2 * hx2

/-- **The mass of a piece**: `∫_r λ_a λ_b = |det β| Area(p)/12 · Σ_{j,l} (1 + δ_jl) β_{ja} β_{lb}`,
the contribution of the fan triangle `r` to `(P₀ᵀ M_S P₀)_{ab}` as `Transfer::new` assembles it.

[formula: ic-transfer-mass-integral] -/
theorem integral_piece_mass {p r : Fin 3 → ℝ × ℝ} (hp : D p ≠ 0) (hr : D r ≠ 0) (a b : Fin 3) :
    ∫ x in triangle r, bary p a x * bary p b x =
      |bdet (pieceRows p r)| * area p / 12 *
        ∑ j : Fin 3, ∑ l : Fin 3, (if j = l then 2 else 1) * pieceRows p r j a * pieceRows p r l b := by
  have e : ∀ x, bary p a x * bary p b x =
      ∑ j, ∑ l, bary p a (r j) * bary p b (r l) * (bary r j x * bary r l x) := by
    intro x
    rw [bary_eq_interp p hr a x, bary_eq_interp p hr b x, interp, interp, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun l _ => by ring
  simp_rw [e]
  have hij : ∀ j l, Integrable (fun x => bary p a (r j) * bary p b (r l) * (bary r j x * bary r l x))
      (volume.restrict (triangle r)) := fun j l =>
    integrableOn_triangle r (continuous_const.mul ((continuous_bary r j).mul (continuous_bary r l)))
  rw [integral_finsetSum _ fun j _ => integrable_finsetSum _ fun l _ => hij j l]
  simp_rw [integral_finsetSum _ fun l _ => hij _ l, integral_const_mul, integral_bary_mul_bary hr,
    area_piece hp r]
  simp only [Fin.sum_univ_three, pieceRows, Fin.isValue, Fin.reduceEq, ite_true, ite_false]
  ring

/-- **The transfer mass matrix is the Galerkin mass matrix of the input triangle.** If the pieces
`r_k` (`k ∈ s`) of the common subdivision tile the input triangle `p`, the per-piece assembly
of `(P₀ᵀ M_S P₀)_{ab}` sums to `∫_p λ_a λ_b = Area(p)(1 + δ_ab)/12`, the entry of `M₀`.

[formula: ic-transfer-mass-integral] -/
theorem transfer_mass_tiling {ι : Type*} {s : Finset ι} {r : ι → Fin 3 → ℝ × ℝ}
    {p : Fin 3 → ℝ × ℝ} (hp : D p ≠ 0) (h : IsTiling s r (triangle p)) (a b : Fin 3) :
    ∑ k ∈ s, |bdet (pieceRows p (r k))| * area p / 12 *
        ∑ j : Fin 3, ∑ l : Fin 3,
          (if j = l then 2 else 1) * pieceRows p (r k) j a * pieceRows p (r k) l b =
      area p * (if a = b then 2 else 1) / 12 := by
  have e := h.integral_eq_sum_of_continuous (f := fun x => bary p a x * bary p b x)
    ((continuous_bary p a).mul (continuous_bary p b))
  rw [Finset.sum_congr rfl fun k hk => (integral_piece_mass hp (h.nondeg k hk) a b).symm, ← e,
    integral_bary_mul_bary hp]
  split_ifs <;> ring

/-- **With pieces of the orientation of `p`** (`det β > 0`, as for the faces of the overlay), the
valuation `pieceMass` of `FEM/Transfer.lean` sums to `(1 + δ_ab)/12 = pieceMass` of the
whole triangle (`pieceMass_one`).

[formula: ic-transfer-mass-integral] -/
theorem sum_pieceMass_tiling {ι : Type*} {s : Finset ι} {r : ι → Fin 3 → ℝ × ℝ}
    {p : Fin 3 → ℝ × ℝ} (hp : D p ≠ 0) (h : IsTiling s r (triangle p))
    (hor : ∀ k ∈ s, 0 < D (r k) * D p) (a b : Fin 3) :
    ∑ k ∈ s, pieceMass (pieceRows p (r k)) a b = (if a = b then 2 else 1) / 12 := by
  have hA : 0 < area p := by
    have := abs_D_eq_two_mul_area p
    have := abs_pos.2 hp
    linarith
  have hdet : ∀ k ∈ s, 0 < bdet (pieceRows p (r k)) := fun k hk => by
    have h1 : 0 < bdet (pieceRows p (r k)) * (D p * D p) := by
      rw [← mul_assoc, ← D_eq_bdet_mul hp]; exact hor k hk
    exact (pos_iff_pos_of_mul_pos h1).2 (mul_self_pos.2 hp)
  have e := transfer_mass_tiling hp h a b
  rw [Finset.sum_congr rfl fun k hk => by rw [abs_of_pos (hdet k hk)]] at e
  have e2 : ∑ k ∈ s, pieceMass (pieceRows p (r k)) a b * area p =
      area p * (if a = b then 2 else 1) / 12 := by
    rw [← e]
    refine Finset.sum_congr rfl fun k _ => ?_
    simp only [pieceMass]
    ring
  rw [← Finset.sum_mul] at e2
  field_simp at e2 ⊢
  linarith

end Transfer

end FEM

end DDG
