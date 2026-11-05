import DdgProofs.Delaunay.Termination

/-!
# Delaunay triangulations, layer 3: the flip loop as the Rust code runs it

`crates/ddg-intrinsic/src/delaunay.rs` (`flip_to_delaunay`, `run_flip_queue`) runs the flip
algorithm of `DdgProofs/Delaunay/Termination.lean` with three differences:

1. **A tolerance.** An edge is flipped when it *fails a criterion* (`fails_delaunay`), not when it
   is merely not locally Delaunay: `AngleSum τ` flips when `cot θ_k + cot θ_l < 0` and
   `(θ_k − π) + θ_l > τ` (`angleSumCriterion`), `CotanWeight τ` when
   `½ (cot θ_k + cot θ_l) < −τ` (`cotanCriterion`). Both are a `FlipCriterion`: a property of the
   five diamond lengths that implies `θ_k + θ_l > π` and does not depend on the halfedge of the
   edge it is read from. So **every flip the loop makes is a Delaunay flip**
   (`IsMetric.delaunayFlippable_of_fails`), and in exact arithmetic the refusals
   `BoundaryEdge`, `DegreeOneEndpoint` (self-glued edges are always Delaunay, layer 3 (b)) and
   `NotConvex` (B–S Prop. 11) never happen (`IsMetric.convexAt_of_fails`). `TooThin` (the shape
   limit of SPEC N5) is a floating-point policy outside this model; constrained edges are not
   modelled (the default has none).
2. **A queue in a fixed order.** All edges start on a first-in first-out queue; an edge is taken
   off, tested, flipped if it fails, and then the four outer edges of its diamond are put back
   (the new edge is not). The order is irrelevant to correctness: `LoopStep` takes *any* queued
   halfedge. The loop invariant (`QueueInv`): **every halfedge that is not queued passes the
   test**. It holds initially (all queued), after a test that passes (`queueInv_pass`), and after
   a flip (`queueInv_flip`), because a flip changes the diamond of no edge except the new edge,
   which is strictly Delaunay (layer 2's `flip_edge_isLocallyDelaunay`), and the four outer
   edges, which are queued (`fails_flip_iff_of_far`: an edge none of whose two faces is one of the
   two faces of the flip keeps its five lengths).
3. **A cap** on the number of flips (`FlipLimit`, default `64·|E|`).

## Main results

* `loop_acc`: if the states reachable from the start by Delaunay flips are finitely many
  (`Termination.lean`; the literature's hypothesis `BoundedLengthFinite` gives it), **every run of
  the loop is finite**, whatever the order (lexicographic in the rank of the state and the queue
  size). Without the cap (`FlipLimit::Unlimited`) the Rust loop therefore exits, after at most
  `|reach s₀| − 1` flips (`MetricGluing.flips_length_le`: the flips of a run are a Delaunay flip
  sequence, and no state repeats, `DelaunayFlip.transGen_ne`).
* `loop_result`: a run from `⟨s₀, all halfedges⟩` that stops (the queue is empty, which is the
  only way to stop, `loopStep_of_mem`) ends in a state in which **no interior edge fails the
  criterion** — the final rescan of `flip_to_delaunay` finds `n_not_delaunay = 0` and reports
  `DelaunayStatus::Delaunay` — with the Euler characteristic, the boundary, the total area and
  the cone angles of `s₀`, and the potential `≤` that of `s₀`. For `AngleSum τ` with `τ ≥ 0` this
  is `θ_k + θ_l ≤ π + τ` at every interior edge (`angleSum_not_fails_iff`); for `τ = 0` it is
  exactly `IsDelaunay` (`isDelaunay_of_angleSum_zero`); for `CotanWeight τ` it is
  `½ (cot θ_k + cot θ_l) ≥ −τ` (by definition).
* A run stopped by the cap (`DelaunayStatus::FlipLimit`) is a prefix of a run: its state is
  reachable from `s₀` by Delaunay flips (`loop_reach`), so it has all the invariants and a lower
  potential, but the theorem says nothing about its Delaunay-ness.

Everything here is exact arithmetic; the floating-point loop is not verified.
-/

namespace DdgProofs.Delaunay

open Real Relation

/-! ## Flip criteria -/

/-- A test on the five diamond lengths `(a, b, c, d, e)` of an edge (layer 1's notation) that the
flip loop uses to decide whether to flip: it may only fire at an edge that is not locally
Delaunay, and it reads the same from both halfedges of the edge (`(a, e, d, c, b)` is the diamond
read from the twin). -/
structure FlipCriterion where
  /-- The edge fails the test and is flipped. -/
  fails : ℝ → ℝ → ℝ → ℝ → ℝ → Prop
  not_locallyDelaunay : ∀ {a b c d e : ℝ}, IsTriangle a b c → IsTriangle a d e →
    fails a b c d e → ¬ IsLocallyDelaunay a b c d e
  mirror : ∀ a b c d e, fails a e d c b ↔ fails a b c d e

theorem cornerCot_swap (a b c : ℝ) : cornerCot a c b = cornerCot a b c := by
  unfold cornerCot
  rw [triArea_swap]
  ring

/-- The Rust test `DelaunayCriterion::AngleSum { tolerance: τ }` (`fails_delaunay`): the sum of the
opposite cotangents is negative and the excess `(θ_k − π) + θ_l` exceeds `τ`.

[formula: it-delaunay-flip-loop] -/
noncomputable def angleSumCriterion (τ : ℝ) : FlipCriterion where
  fails a b c d e :=
    cornerCot a b c + cornerCot a d e < 0 ∧ τ < (cornerAngle a b c - π) + cornerAngle a d e
  not_locallyDelaunay hk hl h := not_isLocallyDelaunay_iff.2 ((pi_lt_iff_cot_neg hk hl).2 h.1)
  mirror a b c d e := by
    rw [cornerCot_swap a e d, cornerCot_swap a c b, Intrinsic.cornerAngle_swap a e d,
      Intrinsic.cornerAngle_swap a c b]
    constructor <;> rintro ⟨h1, h2⟩ <;> exact ⟨by linarith, by linarith⟩

/-- The Rust test `DelaunayCriterion::CotanWeight { tolerance: τ }`: the cotan weight
`½ (cot θ_k + cot θ_l)` is below `−τ`; `τ ≥ 0` (Rust rejects a negative tolerance).

[formula: it-delaunay-flip-loop] -/
noncomputable def cotanCriterion (τ : ℝ) (hτ : 0 ≤ τ) : FlipCriterion where
  fails a b c d e := 0.5 * (cornerCot a b c + cornerCot a d e) < -τ
  not_locallyDelaunay hk hl h := not_isLocallyDelaunay_iff.2 ((pi_lt_iff_cot_neg hk hl).2
    (by linarith))
  mirror a b c d e := by
    rw [cornerCot_swap a e d, cornerCot_swap a c b, add_comm]

/-- **What passing the angle-sum test means**: for `τ ≥ 0`, an edge passes iff
`θ_k + θ_l ≤ π + τ` (the documented guarantee of `DelaunayCriterion::AngleSum`). -/
theorem angleSum_not_fails_iff {τ a b c d e : ℝ} (hτ : 0 ≤ τ) (hk : IsTriangle a b c)
    (hl : IsTriangle a d e) :
    ¬ (angleSumCriterion τ).fails a b c d e ↔ cornerAngle a b c + cornerAngle a d e ≤ π + τ := by
  change ¬ (_ ∧ _) ↔ _
  constructor
  · intro h
    by_contra h'
    exact h ⟨(pi_lt_iff_cot_neg hk hl).1 (by linarith), by linarith⟩
  · rintro h ⟨-, h2⟩
    linarith

/-! ## The test on a triangulation -/

namespace Gluing

variable {F : Type*} {G : Gluing F} {C : FlipCriterion} {ℓ : Halfedge F → ℝ}

/-- The halfedge `h` **fails** the test (Rust `fails_delaunay`): it is interior and its diamond
fails the criterion. Boundary edges never fail. -/
def Fails (C : FlipCriterion) (G : Gluing F) (ℓ : Halfedge F → ℝ) (h : Halfedge F) : Prop :=
  G.twin h ≠ h ∧ C.fails (ℓ h) (ℓ h.prev) (ℓ h.next) (ℓ (G.twin h).next) (ℓ (G.twin h).prev)

/-- The test is a property of the edge. -/
theorem IsMetric.fails_twin (hm : G.IsMetric ℓ) (h : Halfedge F) :
    G.Fails C ℓ (G.twin h) ↔ G.Fails C ℓ h := by
  unfold Fails
  rw [G.twin_twin, hm.twin, C.mirror]
  exact ⟨fun ⟨h1, h2⟩ => ⟨Ne.symm h1, h2⟩, fun ⟨h1, h2⟩ => ⟨Ne.symm h1, h2⟩⟩

theorem IsMetric.not_locallyDelaunayAt_of_fails (hm : G.IsMetric ℓ) {h : Halfedge F}
    (hfl : G.Fails C ℓ h) : ¬ G.LocallyDelaunayAt ℓ h := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h
  exact C.not_locallyDelaunay hk hl hfl.2

/-- **Every edge that fails the test admits a Delaunay flip**: it is not locally Delaunay, and
its two faces are distinct (a self-glued edge is Delaunay). So the Rust refusals `BoundaryEdge`
and `DegreeOneEndpoint` never apply to it.

[formula: it-delaunay-flip-loop] -/
theorem IsMetric.delaunayFlippable_of_fails (hm : G.IsMetric ℓ) {h : Halfedge F}
    (hfl : G.Fails C ℓ h) : G.DelaunayFlippable ℓ h := by
  refine ⟨?_, hm.not_locallyDelaunayAt_of_fails hfl⟩
  by_contra hnf
  exact hm.not_locallyDelaunayAt_of_fails hfl (hm.locallyDelaunayAt_of_not_flippable hfl.1 hnf)

/-- An edge that fails the test has a strictly convex diamond (B–S Prop. 11): the Rust refusal
`NotConvex` never applies to it in exact arithmetic. -/
theorem IsMetric.convexAt_of_fails (hm : G.IsMetric ℓ) {h : Halfedge F} (hfl : G.Fails C ℓ h) :
    G.ConvexAt ℓ h := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h
  exact isConvexDiamond_of_not_isLocallyDelaunay hk hl (hm.not_locallyDelaunayAt_of_fails hfl)

/-- **Locality of the test under a flip**: an edge none of whose faces is a face of the flipped
diamond keeps its twin and its five lengths, hence its test result. -/
theorem fails_flip_iff_of_far [DecidableEq F] {h₀ x : Halfedge F} (hf : G.Flippable h₀)
    (x0 : x.1 ≠ h₀.1) (x1 : x.1 ≠ (G.twin h₀).1) (t0 : (G.twin x).1 ≠ h₀.1)
    (t1 : (G.twin x).1 ≠ (G.twin h₀).1) :
    (G.flip h₀).Fails C (G.flipLen ℓ h₀) x ↔ G.Fails C ℓ x := by
  have hfix : ∀ y : Halfedge F, y.1 ≠ h₀.1 → y.1 ≠ (G.twin h₀).1 → G.flipPerm h₀ y = y :=
    fun y y0 y1 => flipPerm_of_ne (fun e => y0 (by rw [e])) (fun e => y0 (by rw [e]; rfl))
      (fun e => y1 (by rw [e])) (fun e => y1 (by rw [e]; rfl))
  have htw : (G.flip h₀).twin x = G.twin x := by
    have := flip_twin_apply G h₀ x
    rwa [hfix x x0 x1, hfix _ t0 t1] at this
  unfold Fails
  rw [htw, flipLen_of_face hf ℓ x0 x1, flipLen_of_face hf ℓ (x := x.prev) x0 x1,
    flipLen_of_face hf ℓ (x := x.next) x0 x1, flipLen_of_face hf ℓ (x := (G.twin x).next) t0 t1,
    flipLen_of_face hf ℓ (x := (G.twin x).prev) t0 t1]

end Gluing

/-! ## The queue -/

namespace MetricGluing

open Gluing

variable {F : Type*} [DecidableEq F] {C : FlipCriterion}

/-- The halfedges of the four outer edges of the diamond of `h₀` after its flip (Rust: the edges
of `next h`, `prev h`, `next t`, `prev t` for the new edge `h`, `t`), in the numbering of the
flipped triangulation: `h₀`, `prev h₀`, `twin h₀`, `prev (twin h₀)` and their new twins. -/
def outer (G : Gluing F) (h₀ : Halfedge F) : Finset (Halfedge F) :=
  {h₀, h₀.prev, G.twin h₀, (G.twin h₀).prev, (G.flip h₀).twin h₀, (G.flip h₀).twin h₀.prev,
    (G.flip h₀).twin (G.twin h₀), (G.flip h₀).twin (G.twin h₀).prev}

/-- The queue after the flip of `h₀`: the edge of `h₀` leaves it, the other queued halfedges keep
their places under their new ids (`σ`), and the four outer edges are added. -/
def flipQueue (G : Gluing F) (h₀ : Halfedge F) (Q : Finset (Halfedge F)) : Finset (Halfedge F) :=
  ((Q.erase h₀).erase (G.twin h₀)).image (G.flipPerm h₀) ∪ outer G h₀

/-- **The loop invariant**: every halfedge off the queue passes the test. -/
def QueueInv (C : FlipCriterion) (s : MetricGluing F) (Q : Finset (Halfedge F)) : Prop :=
  ∀ h, h ∉ Q → ¬ s.G.Fails C s.ℓ h

omit [DecidableEq F] in
theorem queueInv_univ (s : MetricGluing F) [Fintype F] : QueueInv C s Finset.univ :=
  fun h hh => absurd (Finset.mem_univ h) hh

/-- A test that passes takes the edge off the queue (both halfedges). -/
theorem queueInv_pass {s : MetricGluing F} {Q : Finset (Halfedge F)} {h : Halfedge F}
    (hQ : QueueInv C s Q) (hp : ¬ s.G.Fails C s.ℓ h) :
    QueueInv C s ((Q.erase h).erase (s.G.twin h)) := by
  intro y hy
  simp only [Finset.mem_erase, not_and] at hy
  by_cases e1 : y = s.G.twin h
  · rw [e1, s.isMetric.fails_twin]; exact hp
  by_cases e2 : y = h
  · rw [e2]; exact hp
  exact hQ y fun hyQ => hy e1 e2 hyQ

/-- The state after the flip of a halfedge that fails the test. -/
noncomputable def flipFailing (s : MetricGluing F) (h₀ : Halfedge F)
    (hfl : s.G.Fails C s.ℓ h₀) : MetricGluing F :=
  s.flipAt h₀ (s.isMetric.delaunayFlippable_of_fails hfl)

/-- **The invariant survives a flip.** The new edge is strictly Delaunay, so it passes; the four
outer edges are queued; every other edge keeps its diamond (`fails_flip_iff_of_far`) and its
place on the queue. -/
theorem queueInv_flip {s : MetricGluing F} {Q : Finset (Halfedge F)} {h₀ : Halfedge F}
    (hQ : QueueInv C s Q) (hfl : s.G.Fails C s.ℓ h₀) :
    QueueInv C (s.flipFailing h₀ hfl) (flipQueue s.G h₀ Q) := by
  obtain ⟨hf, hnd⟩ := s.isMetric.delaunayFlippable_of_fails hfl
  set G := s.G
  set G' := G.flip h₀
  set σ := G.flipPerm h₀
  have hm' : G'.IsMetric (G.flipLen s.ℓ h₀) := s.isMetric.flip_of_not_locallyDelaunay hf hnd
  -- the new edge passes
  have hnew : ¬ G'.Fails C (G.flipLen s.ℓ h₀) h₀.next := by
    intro hfails
    refine hm'.not_locallyDelaunayAt_of_fails hfails ((hm'.locallyDelaunayAt_iff h₀.next).2 ?_)
    exact (s.isMetric.flip_edge_isLocallyDelaunay hf hnd).le
  have hnew' : ¬ G'.Fails C (G.flipLen s.ℓ h₀) (G.twin h₀).next := by
    rw [← flip_twin_next hf, hm'.fails_twin]; exact hnew
  have six : ∀ z : Halfedge F, z.1 = h₀.1 ∨ z.1 = (G.twin h₀).1 →
      z = h₀ ∨ z = h₀.next ∨ z = h₀.prev ∨ z = G.twin h₀ ∨ z = (G.twin h₀).next ∨
        z = (G.twin h₀).prev := by
    rintro z (hz | hz)
    · rcases Halfedge.eq_or_eq_next_or_eq_prev hz with h | h | h <;> simp [h]
    · rcases Halfedge.eq_or_eq_next_or_eq_prev hz with h | h | h <;> simp [h]
  have memO : ∀ z ∈ outer G h₀, z ∈ flipQueue G h₀ Q := fun z hz =>
    Finset.mem_union_right _ hz
  have tw : ∀ z, G'.twin (G'.twin z) = z := G'.twin_twin
  intro y hy hfails
  change ¬ y ∈ flipQueue G h₀ Q at hy
  change G'.Fails C (G.flipLen s.ℓ h₀) y at hfails
  by_cases hy0 : y.1 = h₀.1 ∨ y.1 = (G.twin h₀).1
  · rcases six y hy0 with e | e | e | e | e | e <;> subst e
    · exact hy (memO _ (by simp [outer]))
    · exact hnew hfails
    · exact hy (memO _ (by simp [outer]))
    · exact hy (memO _ (by simp [outer]))
    · exact hnew' hfails
    · exact hy (memO _ (by simp [outer]))
  by_cases hy1 : (G'.twin y).1 = h₀.1 ∨ (G'.twin y).1 = (G.twin h₀).1
  · have hyz : y = G'.twin (G'.twin y) := (tw y).symm
    rcases six _ hy1 with e | e | e | e | e | e <;> rw [e] at hyz
    · exact hy (memO _ (by rw [hyz]; simp [outer, G']))
    · rw [flip_twin_next hf] at hyz; rw [hyz] at hfails; exact hnew' hfails
    · exact hy (memO _ (by rw [hyz]; simp [outer, G']))
    · exact hy (memO _ (by rw [hyz]; simp [outer, G']))
    · rw [flip_twin_twin_next hf] at hyz; rw [hyz] at hfails; exact hnew hfails
    · exact hy (memO _ (by rw [hyz]; simp [outer, G']))
  -- far from the diamond
  push Not at hy0 hy1
  obtain ⟨y0, y1⟩ := hy0
  have hfix : ∀ z : Halfedge F, z.1 ≠ h₀.1 → z.1 ≠ (G.twin h₀).1 → σ z = z :=
    fun z z0 z1 => flipPerm_of_ne (fun e => z0 (by rw [e])) (fun e => z0 (by rw [e]; rfl))
      (fun e => z1 (by rw [e])) (fun e => z1 (by rw [e]; rfl))
  have hσy := hfix y y0 y1
  have htw : G'.twin y = σ (G.twin y) := by
    have := flip_twin_apply G h₀ y
    rwa [hσy] at this
  have ht : (G.twin y).1 ≠ h₀.1 ∧ (G.twin y).1 ≠ (G.twin h₀).1 := by
    by_contra hc
    rw [not_and_or, not_not, not_not] at hc
    rcases six _ hc with e | e | e | e | e | e
    · exact y1 (by rw [G.twin_eq_iff.1 e])
    · rw [e, flipPerm_next hf] at htw; exact hy1.2 (by rw [htw])
    · rw [e, flipPerm_of_ne (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)
        (hf.ne h₀.prev (G.twin h₀)) (hf.ne h₀.prev (G.twin h₀).next)] at htw
      exact hy1.1 (by rw [htw]; rfl)
    · exact y0 (by rw [G.twin_eq_iff.1 e, G.twin_twin])
    · rw [e, flipPerm_twin_next hf] at htw; exact hy1.1 (by rw [htw])
    · rw [e, flipPerm_of_ne (hf.ne' (G.twin h₀).prev h₀) (hf.ne' (G.twin h₀).prev h₀.next)
        (Halfedge.prev_ne _) (Halfedge.prev_ne_next _)] at htw
      exact hy1.2 (by rw [htw]; rfl)
  rw [fails_flip_iff_of_far hf y0 y1 ht.1 ht.2] at hfails
  refine hQ y (fun hyQ => hy (Finset.mem_union_left _ (Finset.mem_image.2 ⟨y, ?_, hσy⟩))) hfails
  simp only [Finset.mem_erase]
  exact ⟨fun e => y1 (by rw [e]), fun e => y0 (by rw [e]), hyQ⟩

/-! ## The loop -/

/-- A configuration of the loop: the current triangulation and the queue. -/
structure Config (F : Type*) where
  state : MetricGluing F
  queue : Finset (Halfedge F)

/-- **One step of the loop** (`run_flip_queue`), for any choice of the queued halfedge (Rust takes
the front of a FIFO queue of edges): test it; if it passes, take its edge off the queue; if it
fails, flip it and update the queue (`flipQueue`). -/
inductive LoopStep (C : FlipCriterion) : Config F → Config F → Prop
  | pass {s : MetricGluing F} {Q : Finset (Halfedge F)} {h : Halfedge F} :
      h ∈ Q → ¬ s.G.Fails C s.ℓ h → LoopStep C ⟨s, Q⟩ ⟨s, (Q.erase h).erase (s.G.twin h)⟩
  | flip {s : MetricGluing F} {Q : Finset (Halfedge F)} {h : Halfedge F}
      (hfl : s.G.Fails C s.ℓ h) :
      h ∈ Q → LoopStep C ⟨s, Q⟩ ⟨s.flipFailing h hfl, flipQueue s.G h Q⟩

/-- The loop can always go on while the queue is not empty: it stops exactly when the queue is
empty. -/
theorem loopStep_of_mem {c : Config F} {h : Halfedge F} (hh : h ∈ c.queue) :
    ∃ c', LoopStep C c c' := by
  obtain ⟨s, Q⟩ := c
  by_cases hfl : s.G.Fails C s.ℓ h
  · exact ⟨_, LoopStep.flip hfl hh⟩
  · exact ⟨_, LoopStep.pass hh hfl⟩

/-- A step of the loop keeps the invariant and either keeps the state or makes one Delaunay
flip. -/
theorem LoopStep.inv {c c' : Config F} (hstep : LoopStep C c c') (hQ : QueueInv C c.state c.queue) :
    QueueInv C c'.state c'.queue ∧
      (c'.state = c.state ∧ c'.queue.card < c.queue.card ∨ DelaunayFlip c.state c'.state) := by
  cases hstep with
  | @pass s Q h hh hp =>
    refine ⟨queueInv_pass hQ hp, Or.inl ⟨rfl, ?_⟩⟩
    exact lt_of_le_of_lt (Finset.card_erase_le) (Finset.card_erase_lt_of_mem hh)
  | @flip s Q h hfl hh =>
    exact ⟨queueInv_flip hQ hfl,
      Or.inr (s.delaunayFlip_flipAt h (s.isMetric.delaunayFlippable_of_fails hfl))⟩

/-- **Every configuration a run reaches** has the invariant, and its triangulation is reachable
from the start by Delaunay flips (so it has the start's `χ`, boundary, area and cone angles, and a
lower potential, `MetricGluing.flip_algorithm`). This holds also for a run stopped by the cap. -/
theorem loop_reach {c₀ c : Config F} (h : ReflTransGen (LoopStep C) c₀ c)
    (h0 : QueueInv C c₀.state c₀.queue) :
    QueueInv C c.state c.queue ∧ ReflTransGen DelaunayFlip c₀.state c.state := by
  induction h with
  | refl => exact ⟨h0, ReflTransGen.refl⟩
  | tail _ hst ih =>
    obtain ⟨hinv, hreach⟩ := ih
    obtain ⟨hinv', hor⟩ := hst.inv hinv
    refine ⟨hinv', ?_⟩
    rcases hor with ⟨e, -⟩ | hd
    · rw [e]; exact hreach
    · exact ReflTransGen.tail hreach hd

variable [Fintype F]

/-- **Every run of the loop is finite**, in any order, if finitely many states are reachable from
the start by Delaunay flips (`MetricGluing.reach_finite_of_boundedLengthFinite` gives this from
the literature's hypothesis). A flip lowers the rank of the state among the reachable states, a
test that passes shortens the queue.

[formula: it-delaunay-flip-loop] -/
theorem loop_acc {s₀ : MetricGluing F} (hfin : (reach s₀).Finite) (c : Config F)
    (hc : c.state ∈ reach s₀) : Acc (Function.swap (LoopStep C)) c := by
  classical
  set R := hfin.toFinset
  let rank : MetricGluing F → ℕ := fun t => (R.filter fun u => u.potential < t.potential).card
  let μ : Config F → ℕ × ℕ := fun c => (rank c.state, c.queue.card)
  have wf : WellFounded (InvImage (Prod.Lex (· < ·) (· < ·)) μ) :=
    InvImage.wf μ (WellFounded.prod_lex wellFounded_lt wellFounded_lt)
  refine wf.induction (C := fun c => c.state ∈ reach s₀ → Acc (Function.swap (LoopStep C)) c)
    c (fun c ih hc => Acc.intro c fun c' hstep => ?_) hc
  have hstep' : LoopStep C c c' := hstep
  have hlt : ∀ t u, DelaunayFlip t u → t ∈ reach s₀ → rank u < rank t := by
    intro t u htu ht
    have hu : u ∈ reach s₀ := ReflTransGen.tail ht htu
    have hlt := htu.potential_lt
    refine Finset.card_lt_card ((Finset.ssubset_iff_of_subset ?_).2 ⟨u, ?_, ?_⟩)
    · intro w
      simp only [Finset.mem_filter]
      exact fun ⟨hw, h⟩ => ⟨hw, h.trans hlt⟩
    · simp only [Finset.mem_filter, R, Set.Finite.mem_toFinset]
      exact ⟨hu, hlt⟩
    · simp
  cases hstep' with
  | @pass s Q h hh hp =>
    refine ih _ (Prod.Lex.right _ ?_) hc
    exact lt_of_le_of_lt (Finset.card_erase_le) (Finset.card_erase_lt_of_mem hh)
  | @flip s Q h hfl hh =>
    have hd := s.delaunayFlip_flipAt h (s.isMetric.delaunayFlippable_of_fails hfl)
    exact ih _ (Prod.Lex.left _ _ (hlt _ _ hd hc)) (ReflTransGen.tail hc hd)

/-- **What a finished run of the loop gives** (the Rust loop exiting without hitting the cap):
started from all halfedges queued, a configuration from which the loop cannot go on has an empty
queue, no interior edge fails the test (`flip_to_delaunay` then reports
`DelaunayStatus::Delaunay`), and its triangulation is reachable from the start by Delaunay flips,
with the same Euler characteristic, number of boundary halfedges, total area and cone angles, and
no larger potential.

[formula: it-delaunay-flip-loop] -/
theorem loop_result {s₀ : MetricGluing F} {c : Config F}
    (h : ReflTransGen (LoopStep C) ⟨s₀, Finset.univ⟩ c) (hstop : ∀ c', ¬ LoopStep C c c') :
    c.queue = ∅ ∧ (∀ x, ¬ c.state.G.Fails C c.state.ℓ x) ∧
      ReflTransGen DelaunayFlip s₀ c.state ∧
      c.state.G.eulerChar = s₀.G.eulerChar ∧ c.state.G.numBoundary = s₀.G.numBoundary ∧
      totalArea c.state.ℓ = totalArea s₀.ℓ ∧ c.state.potential ≤ s₀.potential ∧
      ∃ e : s₀.G.Vertex ≃ c.state.G.Vertex,
        ∀ v, c.state.G.coneAngle c.state.ℓ (e v) = s₀.G.coneAngle s₀.ℓ v := by
  obtain ⟨hinv, hreach⟩ := loop_reach h (queueInv_univ s₀)
  have hempty : c.queue = ∅ := by
    rw [Finset.eq_empty_iff_forall_notMem]
    intro x hx
    obtain ⟨c', hc'⟩ := loopStep_of_mem (C := C) hx
    exact hstop c' hc'
  have hc := convexReach_of_reach hreach
  refine ⟨hempty, fun x => hinv x (by rw [hempty]; simp), hreach, eulerChar_eq_of_reach hc,
    numBoundary_eq_of_reach hc, totalArea_eq_of_reach hc, potential_le_of_reach hreach,
    coneAngle_equiv_of_reach hc⟩

omit [DecidableEq F] [Fintype F] in
/-- With the angle-sum test at tolerance `0`, a finished run ends in a **Delaunay**
triangulation. -/
theorem isDelaunay_of_angleSum_zero {s : MetricGluing F}
    (h : ∀ x, ¬ s.G.Fails (angleSumCriterion 0) s.ℓ x) : s.G.IsDelaunay s.ℓ := by
  intro x hx
  obtain ⟨hk, hl⟩ := s.isMetric.diamond_isTriangle x
  have := (angleSum_not_fails_iff le_rfl hk hl).1 (fun hf => h x ⟨hx, hf⟩)
  change _ ≤ π
  linarith

omit [DecidableEq F] [Fintype F] in
/-- With the angle-sum test at tolerance `τ ≥ 0`, a finished run ends with
`θ_k + θ_l ≤ π + τ` at every interior edge. -/
theorem angleSum_le_of_not_fails {τ : ℝ} (hτ : 0 ≤ τ) {s : MetricGluing F}
    (h : ∀ x, ¬ s.G.Fails (angleSumCriterion τ) s.ℓ x) (x : Halfedge F) (hx : s.G.twin x ≠ x) :
    cornerAngleAt s.ℓ x.prev + cornerAngleAt s.ℓ (s.G.twin x).prev ≤ π + τ := by
  obtain ⟨hk, hl⟩ := s.isMetric.diamond_isTriangle x
  rw [← s.isMetric.oppositeAngles_eq]
  exact (angleSum_not_fails_iff hτ hk hl).1 (fun hf => h x ⟨hx, hf⟩)

end MetricGluing

end DdgProofs.Delaunay
