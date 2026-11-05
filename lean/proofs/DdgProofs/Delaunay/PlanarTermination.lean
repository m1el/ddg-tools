import DdgProofs.Delaunay.PlanarFlip

/-!
# Planar Delaunay triangulations, part 3: termination of Lawson's flip algorithm

Lawson's algorithm flips, as long as there is one, an interior edge `ij` (triangles `{i,j,k}`,
`{i,j,l}`) with `l` strictly inside the circumcircle of `i, j, k` (`DelaunayFlip`). Every such flip
applies (`isTriangulation_of_delaunayFlip`, `DdgProofs/Delaunay/PlanarFlip.lean`). Here:

* **The potential.** `liftPotential T = Σ_{t ∈ T} area(t) · Σ_{v ∈ t} |v|²` (`area` is the
  Lebesgue measure of the triangle); it is three times the volume under the lifted triangulation
  (the piecewise linear interpolant of `|x|²`), the potential of the lifting-map proof
  (Edelsbrunner, Ch. 2; for the intrinsic form, layer 1's `liftEnergy` in
  `DdgProofs/Delaunay/Potential.lean`). A flip changes it by exactly
  `−sda(i,j,k) · inCircle i j k l / (2 |sda(i,j,k)|)` (`liftPotential_flip`, from `liftQuad_flip`),
  so **every Delaunay flip strictly decreases it** (`liftPotential_delaunayFlip_lt`).
* **Finitely many triangulations.** A triangulation of `S` is a set of three-element subsets of
  `S` (`IsTriangulation.subset_triangles`), so there are at most `2^(|S| choose 3)` of them
  (`finite_triangulations`).
* **Termination.** The flip relation on triangulations of `S` is well founded
  (`wellFounded_flipStep`); a flip sequence has at most `2^#(triangles S)` steps
  (`flip_sequence_length_le`), there is no infinite one (`not_infinite_flip_sequence`); a
  triangulation in which no Delaunay flip applies is locally Delaunay
  (`isLocallyDelaunay_iff_forall_not_delaunayFlip`), and from every triangulation of `S` the
  algorithm reaches a locally Delaunay triangulation of `S` (`exists_flips_to_locallyDelaunay`).

The existence of **some** triangulation of `S` is a hypothesis of `exists_flips_to_locallyDelaunay`
(and of the corresponding Delaunay statement in `DdgProofs/Delaunay/PlanarLemma.lean`); it is
not proved here (an incremental or sweep construction).

## References

* C. L. Lawson, *Software for C¹ surface interpolation*, in: Mathematical Software III, Academic
  Press, 1977, 161–194.
* R. Sibson, *Locally equiangular triangulations*, The Computer Journal 21(3), 243–245, 1978.
* H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge University Press, 2001,
  Ch. 2 (flipping terminates: the lifted triangulation goes down at every flip).
* M. de Berg, O. Cheong, M. van Kreveld, M. Overmars, *Computational Geometry: Algorithms and
  Applications*, 3rd ed., Springer, 2008, Ch. 9 (Algorithm LegalTriangulation; termination by the
  angle vector instead of a potential).

The bibliographic data were not checked against the works here.
-/

noncomputable section

namespace DdgProofs

namespace Delaunay

namespace Planar

open MeasureTheory Set

/-! ## The lifted potential -/

/-- The area of a triangle (Lebesgue measure of its convex hull). -/
def area (t : Finset (ℝ × ℝ)) : ℝ := volume.real (hull t)

/-- `area(t) · Σ_{v ∈ t} |v|²`. -/
def liftWeight (t : Finset (ℝ × ℝ)) : ℝ := area t * ∑ v ∈ t, nsq v

/-- **The potential** `Σ_{t ∈ T} area(t) · Σ_{v ∈ t} |v|²`. -/
def liftPotential (T : Finset (Finset (ℝ × ℝ))) : ℝ := ∑ t ∈ T, liftWeight t

/-- The area of a triangle is half the absolute signed double area. -/
theorem area_triple {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) : area {a, b, c} = |sda a b c| / 2 := by
  have e := TriangleIntegral.volume_real_triangle (p := ![a, b, c]) h
  rw [TriangleIntegral.triangle_eq_convexHull] at e
  simpa [area, hull_triple, TriangleIntegral.D] using e

theorem liftWeight_triple {a b c : ℝ × ℝ} (h : sda a b c ≠ 0) :
    liftWeight {a, b, c} = liftTri a b c := by
  have hab := sda_ne_zero_left h
  have hbc := sda_ne_zero_right h
  have hac := sda_ne_zero_outer h
  rw [liftWeight, area_triple h, liftTri, Finset.sum_insert (by simp [hab, hac]),
    Finset.sum_insert (by simp [hbc]), Finset.sum_singleton, add_assoc]

theorem finset_ne_of_mem {A B : Finset (ℝ × ℝ)} {x : ℝ × ℝ} (hA : x ∈ A) (hB : x ∉ B) : A ≠ B := by
  rintro rfl; exact hB hA

/-- A new triangle of a flip is not a triangle of the old triangulation. -/
theorem new_not_mem {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))} (hT : IsTriangulation S T)
    {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T) {n : Finset (ℝ × ℝ)}
    (hn : IsTri n) (hsub : hull n ⊆ hull {i, j, k} ∪ hull {i, j, l}) (hn1 : n ≠ {i, j, k})
    (hn2 : n ≠ {i, j, l}) : n ∉ T := by
  intro hnT
  have := disjoint_interior_of_subset hT h1 h2 hnT hn1 hn2 hsub
  obtain ⟨a, b, c, h, rfl⟩ := hn
  obtain ⟨x, hx⟩ := interior_hull3_nonempty h
  rw [← hull_triple] at hx
  exact Set.disjoint_left.1 this hx hx

/-- **The potential changes under a flip by exactly** `−sda(i,j,k) · inCircle i j k l / (2|sda(i,j,k)|)`
(for a strictly convex quadrilateral). -/
theorem liftPotential_flip {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) {i j k l : ℝ × ℝ} (h1 : {i, j, k} ∈ T) (h2 : {i, j, l} ∈ T)
    (hkl : k ≠ l) (hc : sda k l i * sda k l j < 0) :
    liftPotential (flip T i j k l) =
      liftPotential T - sda i j k * inCircle i j k l / (2 * |sda i j k|) := by
  have hk := hT.opposite h1 h2 hkl
  have hD1 := hT.sda_ne_zero h1
  have hD2 := hT.sda_ne_zero h2
  have hni : sda k l i ≠ 0 := by rintro h0; simp [h0] at hc
  have hnj : sda k l j ≠ 0 := by rintro h0; simp [h0] at hc
  have hn1 : sda i k l ≠ 0 := by
    have : sda i k l = sda k l i := by unfold sda; ring
    rwa [this]
  have hn2 : sda j k l ≠ 0 := by
    have : sda j k l = sda k l j := by unfold sda; ring
    rwa [this]
  have hij := sda_ne_zero_left hD1
  have hjk := sda_ne_zero_right hD1
  have hik := sda_ne_zero_outer hD1
  have hjl := sda_ne_zero_right hD2
  have hil := sda_ne_zero_outer hD2
  have hU := flip_hull_union hk hc
  have hsub1 : hull {i, k, l} ⊆ hull {i, j, k} ∪ hull {i, j, l} := by
    rw [hull_triple, hull_triple, hull_triple, ← hU]; exact Set.subset_union_left
  have hsub2 : hull {j, k, l} ⊆ hull {i, j, k} ∪ hull {i, j, l} := by
    rw [hull_triple, hull_triple, hull_triple, ← hU]; exact Set.subset_union_right
  have lnot : l ∉ ({i, j, k} : Finset (ℝ × ℝ)) := by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    exact ⟨hil.symm, hjl.symm, hkl.symm⟩
  have knot : k ∉ ({i, j, l} : Finset (ℝ × ℝ)) := by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    exact ⟨hik.symm, hjk.symm, hkl⟩
  have inot : i ∉ ({j, k, l} : Finset (ℝ × ℝ)) := by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    exact ⟨hij, hik, hil⟩
  have ne11 : ({i, k, l} : Finset (ℝ × ℝ)) ≠ {i, j, k} := finset_ne_of_mem (by simp) lnot
  have ne12 : ({i, k, l} : Finset (ℝ × ℝ)) ≠ {i, j, l} := finset_ne_of_mem (by simp) knot
  have ne21 : ({j, k, l} : Finset (ℝ × ℝ)) ≠ {i, j, k} := finset_ne_of_mem (by simp) lnot
  have ne22 : ({j, k, l} : Finset (ℝ × ℝ)) ≠ {i, j, l} := finset_ne_of_mem (by simp) knot
  have ne12' : ({i, k, l} : Finset (ℝ × ℝ)) ≠ {j, k, l} := finset_ne_of_mem (by simp) inot
  have h12 := IsTriangulation.triple_ne hkl hD2
  have nm1 := new_not_mem hT h1 h2 ⟨i, k, l, hn1, rfl⟩ hsub1 ne11 ne12
  have nm2 := new_not_mem hT h1 h2 ⟨j, k, l, hn2, rfl⟩ hsub2 ne21 ne22
  set R := (T.erase {i, j, k}).erase {i, j, l}
  have hT' : T = insert {i, j, k} (insert {i, j, l} R) := by
    rw [Finset.insert_erase (Finset.mem_erase.2 ⟨h12.symm, h2⟩), Finset.insert_erase h1]
  have hPT : liftPotential T = liftWeight {i, j, k} + liftWeight {i, j, l} + ∑ t ∈ R, liftWeight t := by
    conv_lhs => rw [liftPotential, hT']
    rw [Finset.sum_insert (by simp [R, h12]), Finset.sum_insert (by simp [R]), add_assoc]
  have hPF : liftPotential (flip T i j k l) =
      liftWeight {i, k, l} + liftWeight {j, k, l} + ∑ t ∈ R, liftWeight t := by
    rw [liftPotential, flip, Finset.sum_insert, Finset.sum_insert, add_assoc]
    · simp only [Finset.mem_erase, not_and]; intro _ _; exact nm2
    · simp only [Finset.mem_insert, Finset.mem_erase, not_or, not_and]
      exact ⟨ne12', fun _ _ => nm1⟩
  have hq := liftQuad_flip hk hc
  rw [hPF, hPT, liftWeight_triple hD1, liftWeight_triple hD2, liftWeight_triple hn1,
    liftWeight_triple hn2]
  linarith

/-- **Every Delaunay flip strictly decreases the potential.**

[formula: planar-flip-potential] -/
theorem liftPotential_delaunayFlip_lt {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) (h : DelaunayFlip T T') : liftPotential T' < liftPotential T := by
  obtain ⟨i, j, k, l, h1, h2, hkl, hin, rfl⟩ := h
  have hc := crossing_of_inCircumcircle (hT.opposite h1 h2 hkl) hin
  rw [liftPotential_flip hT h1 h2 hkl hc]
  have hD := hT.sda_ne_zero h1
  have : 0 < sda i j k * inCircle i j k l / (2 * |sda i j k|) :=
    div_pos hin (by positivity)
  linarith

/-! ## Finitely many triangulations -/

/-- The three-element subsets of `S`: the candidate triangles. -/
def triangles (S : Finset (ℝ × ℝ)) : Finset (Finset (ℝ × ℝ)) := S.powerset.filter (·.card = 3)

theorem IsTri.card {t : Finset (ℝ × ℝ)} (ht : IsTri t) : t.card = 3 := by
  obtain ⟨a, b, c, h, rfl⟩ := ht
  exact Finset.card_eq_three.2 ⟨a, b, c, sda_ne_zero_left h, sda_ne_zero_outer h,
    sda_ne_zero_right h, rfl⟩

theorem IsTriangulation.subset_triangles {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) : T ⊆ triangles S := fun t ht =>
  Finset.mem_filter.2 ⟨Finset.mem_powerset.2 (hT.subset t ht), (hT.tri t ht).card⟩

/-- **There are finitely many triangulations of `S`.** -/
theorem finite_triangulations (S : Finset (ℝ × ℝ)) : {T | IsTriangulation S T}.Finite :=
  (triangles S).powerset.finite_toSet.subset fun _ hT =>
    Finset.mem_coe.2 (Finset.mem_powerset.2 (Set.mem_ofPred_eq ▸ hT).subset_triangles)

/-! ## Termination -/

/-- One step of the flip algorithm on triangulations of `S`: `T'` comes from the triangulation `T`
by a Delaunay flip. -/
def FlipStep (S : Finset (ℝ × ℝ)) (T' T : Finset (Finset (ℝ × ℝ))) : Prop :=
  IsTriangulation S T ∧ DelaunayFlip T T'

open Classical in
/-- The number of sets of candidate triangles with a smaller potential. -/
def rank (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) : ℕ :=
  ((triangles S).powerset.filter fun U => liftPotential U < liftPotential T).card

theorem rank_le (S : Finset (ℝ × ℝ)) (T : Finset (Finset (ℝ × ℝ))) :
    rank S T ≤ 2 ^ (triangles S).card := by
  classical
  unfold rank
  rw [← Finset.card_powerset]
  convert Finset.card_filter_le _ _

theorem rank_lt_of_flipStep {S : Finset (ℝ × ℝ)} {T T' : Finset (Finset (ℝ × ℝ))}
    (h : FlipStep S T' T) : rank S T' < rank S T := by
  have hlt := liftPotential_delaunayFlip_lt h.1 h.2
  have hT' := isTriangulation_of_delaunayFlip h.1 h.2
  unfold rank
  apply Finset.card_lt_card
  rw [Finset.ssubset_iff_of_subset]
  · refine ⟨T', ?_, ?_⟩
    · rw [Finset.mem_filter]; exact ⟨Finset.mem_powerset.2 hT'.subset_triangles, hlt⟩
    · rw [Finset.mem_filter]; exact fun hmem => lt_irrefl _ hmem.2
  · intro U hU
    rw [Finset.mem_filter] at hU ⊢
    exact ⟨hU.1, hU.2.trans hlt⟩

/-- **Flipping terminates**: the flip relation on triangulations of `S` is well founded.

[formula: planar-flip-termination] -/
theorem wellFounded_flipStep (S : Finset (ℝ × ℝ)) : WellFounded (FlipStep S) :=
  Subrelation.wf (fun h => rank_lt_of_flipStep h) (measure (rank S)).wf

/-- **A flip sequence has at most `2^#(triangles S)` steps.** -/
theorem flip_sequence_length_le {S : Finset (ℝ × ℝ)} {f : ℕ → Finset (Finset (ℝ × ℝ))} {n : ℕ}
    (h0 : IsTriangulation S (f 0)) (hstep : ∀ m < n, DelaunayFlip (f m) (f (m + 1))) :
    n ≤ 2 ^ (triangles S).card := by
  have key : ∀ m ≤ n, IsTriangulation S (f m) ∧ rank S (f m) + m ≤ rank S (f 0) := by
    intro m
    induction m with
    | zero => intro _; exact ⟨h0, by simp⟩
    | succ m ih =>
      intro hm
      obtain ⟨hTm, hr⟩ := ih (by omega)
      have hs := hstep m (by omega)
      have hlt := rank_lt_of_flipStep ⟨hTm, hs⟩
      exact ⟨isTriangulation_of_delaunayFlip hTm hs, by omega⟩
  have := (key n le_rfl).2
  have := rank_le S (f 0)
  omega

/-- **There is no infinite sequence of Delaunay flips.** -/
theorem not_infinite_flip_sequence {S : Finset (ℝ × ℝ)} {f : ℕ → Finset (Finset (ℝ × ℝ))}
    (h0 : IsTriangulation S (f 0)) : ¬ ∀ m, DelaunayFlip (f m) (f (m + 1)) := by
  intro h
  have := flip_sequence_length_le (n := 2 ^ (triangles S).card + 1) h0 (fun m _ => h m)
  omega

/-- **Lawson's algorithm reaches a locally Delaunay triangulation**: from every triangulation of
`S` some sequence of Delaunay flips leads to a triangulation of `S` that is locally Delaunay
(every maximal sequence does: it is finite, `not_infinite_flip_sequence`, and ends where no flip
applies, `isLocallyDelaunay_iff_forall_not_delaunayFlip`).

[formula: planar-flip-termination] -/
theorem exists_flips_to_locallyDelaunay {S : Finset (ℝ × ℝ)} {T : Finset (Finset (ℝ × ℝ))}
    (hT : IsTriangulation S T) :
    ∃ T', Relation.ReflTransGen DelaunayFlip T T' ∧ IsTriangulation S T' ∧ IsLocallyDelaunay T' := by
  obtain ⟨U, ⟨hTU, hU⟩, hmin⟩ := (wellFounded_flipStep S).has_min
    {U | Relation.ReflTransGen DelaunayFlip T U ∧ IsTriangulation S U} ⟨T, Relation.ReflTransGen.refl, hT⟩
  refine ⟨U, hTU, hU, isLocallyDelaunay_iff_forall_not_delaunayFlip.2 fun V hV => ?_⟩
  exact hmin V ⟨hTU.tail hV, isTriangulation_of_delaunayFlip hU hV⟩ ⟨hU, hV⟩

end Planar

end Delaunay

end DdgProofs
