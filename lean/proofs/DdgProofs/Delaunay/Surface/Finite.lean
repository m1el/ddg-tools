import DdgProofs.Delaunay.Surface.FlipEdge
import DdgProofs.Delaunay.Surface.Fan
import DdgProofs.Delaunay.FlipQueue

/-!
# Delaunay layer 4: finitely many triangulations with bounded edges — termination

The finiteness hypothesis of layer 3 (`MetricGluing.BoundedLengthFinite`,
`DdgProofs/Delaunay/Termination.lean`) is a theorem: **for every metric triangulation `s` and
every `L`, the triangulations in the flip class of `s` whose edges are all at most `L` long are
finitely many** (`MetricGluing.boundedLengthFinite`; here even without "up to isomorphism",
since the faces are a fixed finite type). The proof: `s` is overlaid on itself (`Overlay.refl`),
a convex flip keeps an overlay (`Overlay.flip`, `FlipEdge.lean`), so every `t` in the flip class
is overlaid on `s`; every edge of `t` is then the developed length of a word of `s` with at most
`crossBound s L` crossings (`Overlay.exists_word_le`, `Fan.lean`), and those words are finitely
many (`Word.finite_length_le`).

Consequences, all unconditional (Indermitte–Liebling–Troyanov–Clémençon 2001, Thm. 1;
Bobenko–Springborn 2007, Prop. 12, Thm. 15; Fisher–Springborn–Bobenko–Schröder 2007):

* `MetricGluing.reach_finite`: finitely many states are reachable by Delaunay flips;
* `MetricGluing.acc`, `wellFounded_delaunayFlip`, `not_exists_infinite_flips'`: **Delaunay
  flipping terminates**, from every state, in every order;
* `MetricGluing.exists_delaunay`: **a Delaunay triangulation with the original vertices exists**
  (reachable by Delaunay flips);
* `MetricGluing.flip_algorithm'`: the flip algorithm (layer 3) with no hypothesis;
* `MetricGluing.loop_acc'`: the Rust flip loop (`FlipQueue.lean`) always stops, in any order.

Surfaces with boundary, self-loops, multiple edges and self-glued faces are included.
-/

noncomputable section

namespace DDG.Delaunay.Surface

open DDG.Delaunay.Planar
open Set Relation
open scoped Classical

variable {F : Type*}

/-! ## The identity overlay -/

/-- The identity isomorphism of a state. -/
theorem MetricGluing.isIso_refl (s : MetricGluing F) : MetricGluing.IsIso s s :=
  ⟨⟨Equiv.refl F, fun _ => 0, fun h => by simp⟩, fun h => by simp [Gluing.Iso.toHalfedge]⟩

theorem Cell.tri_id (s : MetricGluing F) (f : F) : Cell.tri s (f, Rigid.id) = sLayout s.ℓ f := by
  funext k; simp [Cell.tri, Tri.map]

/-- `s` overlaid on itself: every face in its standard layout, with itself as its only cell. -/
def Overlay.refl (s : MetricGluing F) : Overlay s s where
  lay f := sLayout s.ℓ f
  lay_dist f k := sLayout_dist s.isMetric f k
  lay_pos f := sLayout_pos s.isMetric f
  cells f := {(f, Rigid.id)}
  cover := by
    intro f y hy
    refine Set.mem_iUnion₂.2 ⟨(f, Rigid.id), Finset.mem_singleton_self _, ?_⟩
    rwa [Cell.tri_id]
  nonjunk := by
    intro f c hc
    rw [Finset.mem_singleton] at hc
    subst hc
    rw [Cell.tri_id, Set.inter_self]
    exact Tri.int_nonempty (sLayout_pos s.isMetric f)
  germ_int := by
    intro f y hy
    refine Or.inl ⟨(f, Rigid.id), ?_, by rwa [Cell.tri_id]⟩
    ext c
    simp only [mem_germAt, Finset.mem_singleton]
    constructor
    · exact fun h => h.1
    · rintro rfl; exact ⟨rfl, by rw [Cell.tri_id]; exact Tri.int_subset_hull hy⟩
  germ_edge := by
    intro x hx y hy
    have hΔ : layΔ s (fun f => sLayout s.ℓ f) x = unfold s x := rfl
    have hac : Cell.across s (x.1, Rigid.id) x.2 = ((s.G.twin x).1, layΔ s (fun f => sLayout s.ℓ f) x) := by
      simp only [Cell.across, hΔ, Rigid.id_comp]
    refine Or.inr ⟨(x.1, Rigid.id), x.2, hx, ?_, by rwa [Cell.tri_id]⟩
    rw [hac]
    ext c
    simp only [mem_germAt, Finset.mem_union, Finset.mem_singleton, Finset.mem_image,
      Finset.mem_insert, Cell.mv]
    constructor
    · rintro ⟨h | ⟨c₀, rfl, rfl⟩, -⟩
      · exact Or.inl h
      · right; simp [Rigid.comp_id]
    · rintro (rfl | rfl)
      · exact ⟨Or.inl rfl, by rw [Cell.tri_id]; exact Tri.side_subset_hull (sLayout_pos s.isMetric _) _ hy⟩
      · refine ⟨Or.inr ⟨((s.G.twin x).1, Rigid.id), rfl, by simp [Rigid.comp_id]⟩, ?_⟩
        rw [← hac]
        have hs := Cell.across_side (s := s) (x.1, Rigid.id) x.2
        rw [Cell.tri_id] at hs
        exact Tri.side_subset_hull (Cell.tri_pos _) _ (by rw [hs]; exact hy)
  germ_bdry := by
    intro x hx y hy
    refine ⟨(x.1, Rigid.id), x.2, hx, ?_, by rwa [Cell.tri_id]⟩
    ext c
    simp only [mem_germAt, Finset.mem_singleton]
    constructor
    · exact fun h => h.1
    · rintro rfl
      exact ⟨rfl, by rw [Cell.tri_id]; exact Tri.side_subset_hull (sLayout_pos s.isMetric _) _ hy⟩
  corner := by
    intro f k c hc _
    rw [Finset.mem_singleton] at hc
    subst hc
    exact ⟨k, by rw [Cell.tri_id]⟩

/-- Every state of the flip class of `s` is overlaid on `s`. -/
theorem exists_overlay [DecidableEq F] (s : MetricGluing F) {t : MetricGluing F}
    (h : t ∈ MetricGluing.flipClass s) : Nonempty (Overlay s t) := by
  induction h with
  | refl => exact ⟨Overlay.refl s⟩
  | @tail b c _ hst ih =>
    have hst' : @MetricGluing.ConvexFlip F (fun a b => Classical.propDecidable (a = b)) b c := by
      convert hst
    exact ih.elim fun O => O.flip hst'

/-! ## Finiteness -/

theorem crossBound_mono [Fintype F] [Nonempty F] (s : MetricGluing F) {L L' : ℝ} (h : L ≤ L') :
    crossBound s L ≤ crossBound s L' := by
  unfold crossBound
  have : 2 * L / switchCost s ≤ 2 * L' / switchCost s :=
    div_le_div_of_nonneg_right (by linarith) (switchCost_pos s).le
  have := Nat.ceil_mono this
  exact Nat.add_le_add_left (Nat.mul_le_mul_left _ (Nat.add_le_add_left this 1)) 1

/-- **Every edge length of a state in the flip class with edges `≤ L` lies in a finite set of
reals depending only on `s` and `L`.** -/
theorem length_mem [Fintype F] [DecidableEq F] [Nonempty F] (s : MetricGluing F) (L : ℝ)
    {t : MetricGluing F} (ht : t ∈ MetricGluing.flipClass s) (hL : ∀ h, t.ℓ h ≤ L)
    (x : Halfedge F) :
    t.ℓ x ∈ (Word.devLength s) '' {w : Word F | w.numCross ≤ crossBound s L} := by
  obtain ⟨O⟩ := exists_overlay s ht
  obtain ⟨w, hw, hwl⟩ := O.exists_word_le x
  exact ⟨w, hw.trans (crossBound_mono s (hL x)), hwl⟩

/-- The states whose lengths lie in a finite set are finitely many. -/
theorem finite_states [Fintype F] (Λ : Set ℝ) (hΛ : Λ.Finite) :
    {t : MetricGluing F | ∀ h, t.ℓ h ∈ Λ}.Finite := by
  let κ : MetricGluing F → (Halfedge F → Halfedge F) × (Halfedge F → ℝ) := fun t => (t.G.twin, t.ℓ)
  have hκ : Function.Injective κ := fun t u e =>
    MetricGluing.ext' (congrArg Prod.fst e) (congrArg Prod.snd e)
  have hsub : {t : MetricGluing F | ∀ h, t.ℓ h ∈ Λ} ⊆
      κ ⁻¹' (Set.univ ×ˢ Set.pi Set.univ fun _ => Λ) := by
    intro t ht
    simp only [Set.mem_preimage, Set.mem_prod, Set.mem_univ, true_and, Set.mem_pi, κ]
    exact fun h _ => ht h
  refine Set.Finite.subset ?_ hsub
  exact Set.Finite.preimage hκ.injOn (Set.finite_univ.prod (Set.Finite.pi fun _ => hΛ))

/-- **The finiteness hypothesis of layer 3 holds for every state**: for every `L`, the states of
the flip class with all edges `≤ L` are finitely many (up to isomorphism, and in fact outright).

[formula: it-delaunay-flip-termination] -/
theorem _root_.DDG.Delaunay.MetricGluing.boundedLengthFinite [Fintype F] [DecidableEq F]
    (s : MetricGluing F) : MetricGluing.BoundedLengthFinite s := by
  intro L
  cases isEmpty_or_nonempty F with
  | inl hF =>
    -- no faces: a single state
    refine ⟨{s}, fun t _ _ => ⟨s, Finset.mem_singleton_self s, ?_⟩⟩
    have : s = t := MetricGluing.ext' (funext fun h => (hF.elim h.1)) (funext fun h => (hF.elim h.1))
    rw [← this]
    exact MetricGluing.isIso_refl s
  | inr hF =>
    have hfin := finite_states (F := F) _
      ((Word.finite_length_le (F := F) (crossBound s L)).image (Word.devLength s))
    refine ⟨hfin.toFinset, fun t ht hL => ⟨t, ?_, MetricGluing.isIso_refl t⟩⟩
    rw [Set.Finite.mem_toFinset]
    exact fun x => length_mem s L ht hL x

end DDG.Delaunay.Surface

namespace DDG.Delaunay.MetricGluing

open Relation

variable {F : Type*} [Fintype F] [DecidableEq F]

/-- **Finitely many states are reachable by Delaunay flips.** -/
theorem reach_finite (s : MetricGluing F) : (reach s).Finite :=
  reach_finite_of_boundedLengthFinite (boundedLengthFinite s)

/-- **Delaunay flipping terminates** from every state: `s` is accessible for the converse of
the flip relation.

[formula: it-delaunay-flip-termination] -/
theorem acc (s : MetricGluing F) : Acc (Function.swap DelaunayFlip) s :=
  acc_of_finite (reach_finite s)

/-- The converse of the Delaunay flip relation is well founded. -/
theorem wellFounded_delaunayFlip : WellFounded (Function.swap (DelaunayFlip (F := F))) :=
  wellFounded_of_finite reach_finite

/-- **There is no infinite sequence of Delaunay flips.**

[formula: it-delaunay-flip-termination] -/
theorem not_exists_infinite_flips' (s : MetricGluing F) :
    ¬ ∃ seq : ℕ → MetricGluing F, seq 0 = s ∧ ∀ n, DelaunayFlip (seq n) (seq (n + 1)) :=
  not_exists_infinite_flips (acc s)

/-- **An intrinsic Delaunay triangulation with the original vertices exists**: one is reachable
from every state by Delaunay flips (Indermitte et al. 2001, Thm. 1; Bobenko–Springborn 2007,
Thm. 15; surfaces with boundary and all degenerate gluings included).

[formula: it-delaunay-existence-uniqueness] -/
theorem exists_delaunay (s : MetricGluing F) :
    ∃ t, ReflTransGen DelaunayFlip s t ∧ t.G.IsDelaunay t.ℓ :=
  exists_delaunay_of_acc (acc s)

/-- **The flip algorithm, unconditional** (layer 3's `flip_algorithm` without its hypothesis).

[formula: it-delaunay-flip-algorithm] -/
theorem flip_algorithm' (s : MetricGluing F) :
    Acc (Function.swap DelaunayFlip) s ∧
      (∃ t, ReflTransGen DelaunayFlip s t ∧ t.G.IsDelaunay t.ℓ) ∧
      ∀ t, ReflTransGen DelaunayFlip s t → (∀ u, ¬ DelaunayFlip t u) →
        t.G.IsDelaunay t.ℓ ∧ t.G.eulerChar = s.G.eulerChar ∧
          t.G.numBoundary = s.G.numBoundary ∧ Gluing.totalArea t.ℓ = Gluing.totalArea s.ℓ ∧
          ∃ e : s.G.Vertex ≃ t.G.Vertex, ∀ v, t.G.coneAngle t.ℓ (e v) = s.G.coneAngle s.ℓ v :=
  flip_algorithm (reach_finite s)

/-- A bound on the number of Delaunay flips from `s`: at most `|reach s| − 1`. -/
theorem flips_length_le' (s : MetricGluing F) (seq : ℕ → MetricGluing F) (n : ℕ)
    (h0 : seq 0 = s) (hs : ∀ i < n, DelaunayFlip (seq i) (seq (i + 1))) :
    n + 1 ≤ (reach_finite s).toFinset.card :=
  flips_length_le (reach_finite s) seq n h0 hs

/-- **The Rust flip loop always stops**, in any order (`FlipQueue.lean`'s `loop_acc` without its
hypothesis).

[formula: it-delaunay-flip-loop] -/
theorem loop_acc' {C : FlipCriterion} (s₀ : MetricGluing F) (c : Config F)
    (hc : c.state ∈ reach s₀) : Acc (Function.swap (LoopStep C)) c :=
  loop_acc (reach_finite s₀) c hc

end DDG.Delaunay.MetricGluing
