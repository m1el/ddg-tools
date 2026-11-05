import DdgProofs.Delaunay.Metric
import DdgProofs.Delaunay.Potential
import Mathlib.Data.Set.Finite.Basic
import Mathlib.Data.Fintype.Perm
import Mathlib.Order.WellFounded

/-!
# Delaunay triangulations, layer 3: the flip algorithm on a triangulated surface

The edge-flip algorithm for intrinsic Delaunay triangulations of a piecewise flat surface
(Indermitte, Liebling, Troyanov, Clémençon, *Voronoi diagrams on piecewise flat surfaces and an
application to biological growth*, TCS 263, 2001, Sec. 4, Thm. 1; Bobenko, Springborn, *A discrete
Laplace–Beltrami operator for simplicial surfaces*, DCG 38, 2007, Sec. 2–3, Prop. 10–12, Thm. 15;
Fisher, Springborn, Bobenko, Schröder, *An algorithm for the construction of intrinsic Delaunay
triangulations with applications to digital geometry processing*, Computing 81, 2007, Sec. 2):
*while some interior edge is not locally Delaunay, flip it.* The triangulations are the gluing
data of `DdgProofs/Delaunay/Complex.lean` (self-loops, multiple edges and faces glued to
themselves allowed), with an edge-length metric (`DdgProofs/Delaunay/Metric.lean`).

## Main results

* **Delaunay** (`Gluing.IsDelaunay`): every *interior* edge (`twin h ≠ h`) is locally Delaunay,
  `θ_k + θ_l ≤ π` (B–S Def. 8, Prop. 10; boundary edges impose nothing, as in B–S).
* **Self-glued edges** (`IsMetric.locallyDelaunayAt_of_same_face`): an interior edge whose two
  halfedges lie in the *same* face is always locally Delaunay, strictly — its two opposite angles
  are two distinct corners of one triangle. So an edge on which the flip is undefined
  (`¬ Flippable`, Rust `FlipRefusal::DegreeOneEndpoint`) never needs flipping, and
  `IsMetric.isDelaunay_iff`: **a metric triangulation is Delaunay iff no edge admits a Delaunay
  flip** (an edge that is flippable and not locally Delaunay).
* **The global potential** `P = Σ_faces A(a² + b² + c²)` (`liftPotential`; the lifting energy of
  layer 1, `it-delaunay-potential`) changes only in the two faces of a flip
  (`IsMetric.facePotential_flip`), so it **decreases strictly with every Delaunay flip**
  (`IsMetric.liftPotential_flip_lt`) and is unchanged by the flip of a tie
  (`IsMetric.liftPotential_flip_eq_of_tie`). The same holds for Musin's harmonic index
  `H = Σ (a² + b² + c²)/A` (B–S Prop. 12, `IsMetric.harmonicPotential_flip_lt`), and every
  flip of a strictly convex diamond keeps the total area (`IsMetric.totalArea_flip`).
* **Gauss–Bonnet with boundary** (`IsMetric.gauss_bonnet`): `Σ_v (2π − π b_v − Θ_v) = 2π χ`,
  `b_v` the number of boundary halfedges at `v` (layer 2 had the closed case).
* **States and the flip relation** (`MetricGluing`, `MetricGluing.DelaunayFlip`): a state is
  gluing data with a metric on a fixed face type `F` (a flip keeps the faces). `DelaunayFlip s t`:
  `t` is `s` with one Delaunay flip. Along any sequence of flips of strictly convex diamonds
  (`ConvexFlip`, which every Delaunay flip is) the Euler characteristic, the number of boundary
  halfedges, the boundary halfedges with their lengths, the total area and the cone angle of every
  vertex are unchanged (`ConvexReach.eulerChar_eq`, `numBoundary_eq`, `boundary_equiv`,
  `totalArea_eq`, `coneAngle_equiv`).
* **Termination, conditional** (`MetricGluing.acc_of_finite`): if the states reachable from `s` by
  Delaunay flips are finitely many, then every sequence of Delaunay flips from `s` is finite
  (`Acc`; `not_exists_infinite_flips`), it has at most `|reach s| − 1` flips
  (`flips_length_le`), no state repeats (`DelaunayFlip.transGen_ne`, unconditionally), a
  Delaunay state is reachable (`exists_delaunay_of_acc`), and a state from which no Delaunay flip
  is possible — the end of every maximal sequence — is Delaunay
  (`isDelaunay_iff_forall_not_delaunayFlip`), with the invariants above.
* **The finiteness hypothesis in the form of the literature** (`BoundedLengthFinite`): for every
  `L`, the triangulations in the flip class of `s` (reachable by flips of strictly convex
  diamonds, i.e. geodesic triangulations of the same surface with the same vertices) whose edges
  are all at most `L` long are finitely many *up to isomorphism* (`MetricGluing.Iso`: a
  relabelling of the faces and corners, carrying the lengths). This is the statement that "there
  are only finitely many geodesic triangulations with bounded edge lengths" on a compact
  piecewise flat surface (B–S, proof of Thm. 15 / Prop. 12; Indermitte et al., proof of Thm. 1,
  in the form "the potential takes values in a discrete set"), which rests on the geometry of
  geodesics on the cone surface and was cited in layer 3; **layer 4 proves it for every
  `MetricGluing`** (`Surface/Finite.lean`, `MetricGluing.boundedLengthFinite`, which also gives
  the unconditional `acc`, `exists_delaunay`, `flip_algorithm'`). It implies the finiteness of the reachable set (`reach_finite_of_boundedLengthFinite`):
  along Delaunay flips the harmonic index does not increase and the area is constant, and an edge
  of length `ℓ` in a face of area `A` has `ℓ² ≤ H_face · A ≤ H · A_total`.

**Why finiteness is needed.** A strictly decreasing real potential does not by itself stop a
sequence (`1/n` decreases forever). Both literature proofs use the finiteness of the set of
triangulations below a bound; on a surface there is no bound on the number of flips in terms of
the number of vertices (Indermitte et al., after Thm. 1).

The algorithm as the Rust code runs it (a queue of edges, a tolerance, a cap) is in
`DdgProofs/Delaunay/FlipQueue.lean`.
-/

namespace DDG.Delaunay

open Real

/-! ## Symmetric functions of a triangle -/

/-- A function of the three side lengths of a triangle that does not depend on their order. -/
structure IsSymm3 (φ : ℝ → ℝ → ℝ → ℝ) : Prop where
  rotate : ∀ a b c, φ b c a = φ a b c
  swap : ∀ a b c, φ a c b = φ a b c

theorem isSymm3_triArea : IsSymm3 triArea := ⟨triArea_rotate, triArea_swap⟩

theorem isSymm3_liftEnergy : IsSymm3 liftEnergy :=
  ⟨fun a b c => by unfold liftEnergy; rw [triArea_rotate]; ring,
    fun a b c => by unfold liftEnergy; rw [triArea_swap]; ring⟩

theorem isSymm3_harmonicIndex : IsSymm3 harmonicIndex :=
  ⟨fun a b c => by unfold harmonicIndex; rw [triArea_rotate]; ring_nf,
    fun a b c => by unfold harmonicIndex; rw [triArea_swap]; ring_nf⟩

section
variable {F : Type*} {φ : ℝ → ℝ → ℝ → ℝ}

/-- A symmetric function of a face can be read from any of its halfedges. -/
theorem IsSymm3.face_eq (hφ : IsSymm3 φ) (ℓ : Halfedge F → ℝ) (h : Halfedge F) :
    φ (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) = φ (ℓ h) (ℓ h.next) (ℓ h.prev) := by
  obtain ⟨f, k⟩ := h
  fin_cases k
  · rfl
  · exact (hφ.rotate _ _ _).symm
  · exact (hφ.rotate _ _ _).trans (hφ.rotate _ _ _) |>.symm

theorem IsSymm3.perm_bca (hφ : IsSymm3 φ) (a b c : ℝ) : φ b c a = φ a b c := hφ.rotate a b c
theorem IsSymm3.perm_cab (hφ : IsSymm3 φ) (a b c : ℝ) : φ c a b = φ a b c :=
  (hφ.rotate c a b).symm
theorem IsSymm3.perm_acb (hφ : IsSymm3 φ) (a b c : ℝ) : φ a c b = φ a b c := hφ.swap a b c
theorem IsSymm3.perm_cba (hφ : IsSymm3 φ) (a b c : ℝ) : φ c b a = φ a b c := by
  rw [hφ.swap c a b, hφ.perm_cab]
theorem IsSymm3.perm_bac (hφ : IsSymm3 φ) (a b c : ℝ) : φ b a c = φ a b c := by
  rw [hφ.rotate c b a, hφ.perm_cba]

end

namespace Gluing

variable {F : Type*} {G : Gluing F} {φ : ℝ → ℝ → ℝ → ℝ}

/-! ## Face potentials -/

/-- The sum over the faces of a symmetric function `φ` of the three side lengths. It depends on
the metric only (and on `F`), not on the gluing. -/
noncomputable def facePotential [Fintype F] (φ : ℝ → ℝ → ℝ → ℝ) (ℓ : Halfedge F → ℝ) : ℝ :=
  ∑ f : F, φ (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2))

/-- **The global lifting potential** `P = Σ_faces A (a² + b² + c²)` (layer 1's `liftEnergy`). -/
noncomputable abbrev liftPotential [Fintype F] (ℓ : Halfedge F → ℝ) : ℝ :=
  facePotential liftEnergy ℓ

/-- **Musin's harmonic index** of the triangulation, `H = Σ_faces (a² + b² + c²)/A`. -/
noncomputable abbrev harmonicPotential [Fintype F] (ℓ : Halfedge F → ℝ) : ℝ :=
  facePotential harmonicIndex ℓ

/-- The total area `Σ_faces A`. -/
noncomputable abbrev totalArea [Fintype F] (ℓ : Halfedge F → ℝ) : ℝ := facePotential triArea ℓ

/-! ## Delaunay triangulations -/

section Delaunay

variable {ℓ : Halfedge F → ℝ}

/-- The three corner angles of a face add up to `π`, read from any of its halfedges. -/
theorem IsMetric.cornerAngleAt_add_next_add_prev (hm : G.IsMetric ℓ) (h : Halfedge F) :
    cornerAngleAt ℓ h + cornerAngleAt ℓ h.next + cornerAngleAt ℓ h.prev = π := by
  have t := hm.tri_at h
  have := cornerAngle_add_cornerAngle_add_cornerAngle t
  simp only [cornerAngleAt, Halfedge.next_next, Halfedge.prev_next, Halfedge.next_prev,
    Halfedge.prev_prev]
  rw [Intrinsic.cornerAngle_swap (ℓ h.next) (ℓ h.prev) (ℓ h),
    Intrinsic.cornerAngle_swap (ℓ h.prev) (ℓ h) (ℓ h.next),
    Intrinsic.cornerAngle_swap (ℓ h) (ℓ h.next) (ℓ h.prev)]
  linarith

/-- Every corner angle of a metric is positive. -/
theorem IsMetric.cornerAngleAt_pos (hm : G.IsMetric ℓ) (h : Halfedge F) :
    0 < cornerAngleAt ℓ h :=
  cornerAngle_pos (hm.tri_at h).rotate.swap

/-- The two angles opposite the edge of `h₀`, in the five diamond lengths of layer 1. -/
theorem IsMetric.oppositeAngles_eq (hm : G.IsMetric ℓ) (h₀ : Halfedge F) :
    cornerAngle (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) +
        cornerAngle (ℓ h₀) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev) =
      cornerAngleAt ℓ h₀.prev + cornerAngleAt ℓ (G.twin h₀).prev := by
  simp only [cornerAngleAt, Halfedge.next_prev, Halfedge.prev_prev, hm.twin]
  rw [Intrinsic.cornerAngle_swap (ℓ h₀) (ℓ (G.twin h₀).prev)]

/-- Being locally Delaunay is a property of the edge: it does not depend on which of its two
halfedges is used. -/
theorem IsMetric.locallyDelaunayAt_twin (hm : G.IsMetric ℓ) (h : Halfedge F) :
    G.LocallyDelaunayAt ℓ (G.twin h) ↔ G.LocallyDelaunayAt ℓ h := by
  rw [hm.locallyDelaunayAt_iff, hm.locallyDelaunayAt_iff, G.twin_twin, add_comm]

/-- **An interior edge with the same face on both sides is locally Delaunay, strictly.** If
`twin h ≠ h` lies in the face of `h`, the two corners opposite the edge, `prev h` and
`prev (twin h)`, are two distinct corners of that one triangle (`twin h = next h` gives the corners
`prev h` and `h`, `twin h = prev h` the corners `prev h` and `next h`), so they sum to `π` minus the
third corner angle. This is the edge with an end point of degree one (Rust
`FlipRefusal::DegreeOneEndpoint`; Sharp–Soliman–Crane 2019, Fig. 6), which the flip leaves alone.

[formula: it-delaunay-self-glued] -/
theorem IsMetric.locallyDelaunayAt_of_same_face (hm : G.IsMetric ℓ) {h : Halfedge F}
    (hb : G.twin h ≠ h) (hs : (G.twin h).1 = h.1) :
    cornerAngleAt ℓ h.prev + cornerAngleAt ℓ (G.twin h).prev < π := by
  have hsum := hm.cornerAngleAt_add_next_add_prev h
  rcases Halfedge.eq_or_eq_next_or_eq_prev hs with e | e | e
  · exact absurd e hb
  · rw [e, Halfedge.prev_next]
    have := hm.cornerAngleAt_pos h.next
    linarith
  · rw [e, Halfedge.prev_prev]
    have := hm.cornerAngleAt_pos h
    linarith

/-- A triangulation is **Delaunay** when every interior edge is locally Delaunay
(Bobenko–Springborn 2007, Def. 8 and Prop. 10; Indermitte et al. 2001, Def. 5). Boundary edges
impose no condition. -/
def IsDelaunay (G : Gluing F) (ℓ : Halfedge F → ℝ) : Prop :=
  ∀ h, G.twin h ≠ h → G.LocallyDelaunayAt ℓ h

/-- An edge **admits a Delaunay flip**: the flip is defined (two distinct faces) and the edge is
not locally Delaunay. -/
def DelaunayFlippable (G : Gluing F) (ℓ : Halfedge F → ℝ) (h : Halfedge F) : Prop :=
  G.Flippable h ∧ ¬ G.LocallyDelaunayAt ℓ h

/-- An interior edge on which the flip is not defined is locally Delaunay.

[formula: it-delaunay-self-glued] -/
theorem IsMetric.locallyDelaunayAt_of_not_flippable (hm : G.IsMetric ℓ) {h : Halfedge F}
    (hb : G.twin h ≠ h) (hnf : ¬ G.Flippable h) : G.LocallyDelaunayAt ℓ h :=
  (hm.locallyDelaunayAt_iff h).2 (hm.locallyDelaunayAt_of_same_face hb (not_not.1 hnf)).le

/-- **No edge admits a Delaunay flip iff the triangulation is Delaunay**: an interior edge that is
not locally Delaunay always has two distinct faces (`locallyDelaunayAt_of_not_flippable`), so the
flip algorithm can only stop at a Delaunay triangulation.

[formula: it-delaunay-self-glued] -/
theorem IsMetric.isDelaunay_iff (hm : G.IsMetric ℓ) :
    G.IsDelaunay ℓ ↔ ∀ h, ¬ G.DelaunayFlippable ℓ h := by
  constructor
  · rintro hD h ⟨hf, hnd⟩
    exact hnd (hD h hf.twin_ne)
  · intro hno h hb
    by_cases hf : G.Flippable h
    · by_contra hnd
      exact hno h ⟨hf, hnd⟩
    · exact hm.locallyDelaunayAt_of_not_flippable hb hf

end Delaunay

/-! ## The global potential under a flip -/

section PotentialFlip

variable [Fintype F] [DecidableEq F] {ℓ : Halfedge F → ℝ} {h₀ : Halfedge F}

/-- **A flip changes a face potential only in its two faces**: with the diamond lengths
`a = ℓ h₀`, `b = ℓ (prev h₀)`, `c = ℓ (next h₀)`, `d = ℓ (next (twin h₀))`,
`e = ℓ (prev (twin h₀))` and the new length `f`,
`Φ(ℓ') + φ(a, b, c) + φ(a, d, e) = Φ(ℓ) + φ(f, b, d) + φ(f, e, c)`. -/
theorem IsMetric.facePotential_flip (hφ : IsSymm3 φ) (hm : G.IsMetric ℓ) (hf : G.Flippable h₀) :
    facePotential φ (G.flipLen ℓ h₀) +
        (φ (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next) + φ (ℓ h₀) (ℓ (G.twin h₀).next) (ℓ (G.twin h₀).prev)) =
      facePotential φ ℓ +
        (φ (G.newLength ℓ h₀) (ℓ h₀.prev) (ℓ (G.twin h₀).next) +
          φ (G.newLength ℓ h₀) (ℓ (G.twin h₀).prev) (ℓ h₀.next)) := by
  have hne : (G.twin h₀).1 ≠ h₀.1 := hf
  set R := (Finset.univ.erase h₀.1).erase (G.twin h₀).1
  have key : ∀ ψ : Halfedge F → ℝ, facePotential φ ψ = φ (ψ h₀) (ψ h₀.next) (ψ h₀.prev) +
      (φ (ψ (G.twin h₀)) (ψ (G.twin h₀).next) (ψ (G.twin h₀).prev) +
        ∑ f ∈ R, φ (ψ (f, 0)) (ψ (f, 1)) (ψ (f, 2))) := by
    intro ψ
    rw [← hφ.face_eq ψ h₀, ← hφ.face_eq ψ (G.twin h₀),
      Finset.add_sum_erase (Finset.univ.erase h₀.1) (fun f => φ (ψ (f, 0)) (ψ (f, 1)) (ψ (f, 2)))
        (Finset.mem_erase.2 ⟨hne, Finset.mem_univ _⟩),
      Finset.add_sum_erase Finset.univ (fun f => φ (ψ (f, 0)) (ψ (f, 1)) (ψ (f, 2)))
        (Finset.mem_univ _)]
    rfl
  have rest : ∑ f ∈ R, φ (G.flipLen ℓ h₀ (f, 0)) (G.flipLen ℓ h₀ (f, 1)) (G.flipLen ℓ h₀ (f, 2)) =
      ∑ f ∈ R, φ (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2)) := by
    refine Finset.sum_congr rfl fun f hfR => ?_
    have h1 : f ≠ (G.twin h₀).1 := (Finset.mem_erase.1 hfR).1
    have h2 : f ≠ h₀.1 := (Finset.mem_erase.1 (Finset.mem_erase.1 hfR).2).1
    rw [flipLen_of_face hf ℓ (x := (f, 0)) h2 h1, flipLen_of_face hf ℓ (x := (f, 1)) h2 h1,
      flipLen_of_face hf ℓ (x := (f, 2)) h2 h1]
  rw [key, key, rest, flipLen_self hf, flipLen_next, flipLen_prev hf, flipLen_twin hf,
    flipLen_twin_next, flipLen_twin_prev hf, hm.twin,
    hφ.swap (ℓ h₀) (ℓ h₀.prev) (ℓ h₀.next),
    hφ.perm_cab (G.newLength ℓ h₀) (ℓ h₀.prev) (ℓ (G.twin h₀).next),
    hφ.perm_cab (G.newLength ℓ h₀) (ℓ (G.twin h₀).prev) (ℓ h₀.next)]
  ring

/-- **The global lifting potential decreases strictly with every Delaunay flip** (layer 1's
`liftEnergy_flip_lt` in the two faces of the flip, nothing changes elsewhere).

[formula: it-delaunay-potential] -/
theorem IsMetric.liftPotential_flip_lt (hm : G.IsMetric ℓ) (hf : G.Flippable h₀)
    (hnd : ¬ G.LocallyDelaunayAt ℓ h₀) : liftPotential (G.flipLen ℓ h₀) < liftPotential ℓ := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  have := liftEnergy_flip_lt hk hl hnd
  have e := hm.facePotential_flip isSymm3_liftEnergy hf
  unfold newLength at e
  linarith

/-- At a tie (`θ_k + θ_l = π`) the flip leaves the global lifting potential unchanged. -/
theorem IsMetric.liftPotential_flip_eq_of_tie (hm : G.IsMetric ℓ) (hf : G.Flippable h₀)
    (htie : cornerAngleAt ℓ h₀.prev + cornerAngleAt ℓ (G.twin h₀).prev = π) :
    liftPotential (G.flipLen ℓ h₀) = liftPotential ℓ := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  have := liftEnergy_flip_eq_of_tie hk hl (by rw [hm.oppositeAngles_eq]; exact htie)
  have e := hm.facePotential_flip isSymm3_liftEnergy hf
  unfold newLength at e
  linarith

/-- **Musin's harmonic index decreases strictly with every Delaunay flip** (B–S Prop. 12;
layer 1's `harmonicIndex_flip_lt`).

[formula: it-delaunay-potential] -/
theorem IsMetric.harmonicPotential_flip_lt (hm : G.IsMetric ℓ) (hf : G.Flippable h₀)
    (hnd : ¬ G.LocallyDelaunayAt ℓ h₀) :
    harmonicPotential (G.flipLen ℓ h₀) < harmonicPotential ℓ := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  have := harmonicIndex_flip_lt hk hl hnd
  have e := hm.facePotential_flip isSymm3_harmonicIndex hf
  have s := isSymm3_harmonicIndex.swap (G.newLength ℓ h₀) (ℓ h₀.next) (ℓ (G.twin h₀).prev)
  unfold newLength at e s
  linarith

/-- **A flip of a strictly convex diamond keeps the total area.** -/
theorem IsMetric.totalArea_flip (hm : G.IsMetric ℓ) (hf : G.Flippable h₀)
    (hconv : G.ConvexAt ℓ h₀) : totalArea (G.flipLen ℓ h₀) = totalArea ℓ := by
  obtain ⟨hk, hl⟩ := hm.diamond_isTriangle h₀
  have := flip_area hk hl hconv
  have e := hm.facePotential_flip isSymm3_triArea hf
  unfold newLength at e
  linarith

end PotentialFlip

/-! ## Gauss–Bonnet with boundary -/

/-- The number of boundary halfedges whose tail is the vertex `v` (`1` at a vertex of the boundary
of a surface, `0` at an interior vertex). -/
noncomputable def boundaryCount [Fintype F] [DecidableEq F] (G : Gluing F) (v : G.Vertex) : ℕ :=
  (Finset.univ.filter fun h : Halfedge F => G.twin h = h ∧ G.vertexOf h = v).card

/-- **Discrete Gauss–Bonnet with boundary** on gluing data:
`Σ_v (2π − π b_v − Θ_v) = 2π χ`, with `b_v` the number of boundary halfedges at `v`. At an
interior vertex the term is the angle defect `2π − Θ_v`; at a vertex with one boundary halfedge
it is the boundary turning `π − Θ_v` (the convention of `GaussBonnet.lean`). For a closed surface
this is `gauss_bonnet_closed`.

[formula: gauss-bonnet] -/
theorem IsMetric.gauss_bonnet [Fintype F] [DecidableEq F] {ℓ : Halfedge F → ℝ}
    (hm : G.IsMetric ℓ) :
    ∑ v, (2 * π - π * G.boundaryCount v - G.coneAngle ℓ v) = 2 * π * G.eulerChar := by
  classical
  have hB : ∑ v, (G.boundaryCount v : ℝ) = G.numBoundary := by
    have := Finset.card_eq_sum_card_fiberwise (f := G.vertexOf)
      (s := Finset.univ.filter fun h : Halfedge F => G.twin h = h)
      (t := (Finset.univ : Finset G.Vertex)) (fun _ _ => Finset.mem_univ _)
    rw [numBoundary, this]
    push_cast
    refine Finset.sum_congr rfl fun v _ => ?_
    unfold boundaryCount
    rw [Finset.filter_filter]
  have hE := G.two_mul_numEdges
  have hE' : (2 : ℝ) * G.numEdges = 3 * Fintype.card F + G.numBoundary := by exact_mod_cast hE
  rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, hm.sum_coneAngle,
    ← Finset.mul_sum Finset.univ (fun v => (G.boundaryCount v : ℝ)) π, hB,
    Finset.sum_const, Finset.card_univ, eulerChar, numVertices, Nat.card_eq_fintype_card,
    nsmul_eq_mul]
  push_cast
  linear_combination π * hE'

/-! ## Edge lengths are bounded by the harmonic index and the area -/

/-- **Every edge is bounded by the harmonic index and the total area**: `ℓ(h)² ≤ H · A`, since
`ℓ(h)² ≤ a² + b² + c² = H_face · A_face` for the face of `h`, and every face term of `H` and of
`A` is positive. -/
theorem IsMetric.sq_le_harmonicPotential_mul_totalArea [Fintype F] {ℓ : Halfedge F → ℝ}
    (hm : G.IsMetric ℓ) (h : Halfedge F) : ℓ h ^ 2 ≤ harmonicPotential ℓ * totalArea ℓ := by
  have hpos : ∀ f : F, 0 < triArea (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2)) := fun f => triArea_pos (hm.tri f)
  have hHpos : ∀ f : F, 0 < harmonicIndex (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2)) := fun f => by
    have := (hm.tri f).pos_a
    unfold harmonicIndex
    have := hpos f
    positivity
  have hH : harmonicIndex (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) ≤ harmonicPotential ℓ :=
    Finset.single_le_sum (f := fun f : F => harmonicIndex (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2)))
      (fun f _ => (hHpos f).le) (Finset.mem_univ h.1)
  have hA : triArea (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) ≤ totalArea ℓ :=
    Finset.single_le_sum (f := fun f : F => triArea (ℓ (f, 0)) (ℓ (f, 1)) (ℓ (f, 2)))
      (fun f _ => (hpos f).le) (Finset.mem_univ h.1)
  have hprod : harmonicIndex (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) *
      triArea (ℓ (h.1, 0)) (ℓ (h.1, 1)) (ℓ (h.1, 2)) = ℓ h ^ 2 + ℓ h.next ^ 2 + ℓ h.prev ^ 2 := by
    rw [isSymm3_harmonicIndex.face_eq, isSymm3_triArea.face_eq]
    have := (hpos h.1).ne'
    rw [isSymm3_triArea.face_eq] at this
    unfold harmonicIndex
    field_simp
  have := mul_le_mul hH hA (hpos h.1).le ((hHpos h.1).le.trans hH)
  nlinarith [sq_nonneg (ℓ h.next), sq_nonneg (ℓ h.prev)]

end Gluing

/-! ## States and the flip relation -/

/-- A **triangulated piecewise flat surface**: gluing data on the faces `F` with an edge-length
metric. A flip keeps the face type, so all states of the algorithm live in one type. -/
structure MetricGluing (F : Type*) where
  /-- The gluing data (the combinatorics). -/
  G : Gluing F
  /-- The edge lengths. -/
  ℓ : Halfedge F → ℝ
  isMetric : G.IsMetric ℓ

namespace MetricGluing

open Relation Gluing

variable {F : Type*} {s t : MetricGluing F}

/-- Two states are equal when their twin maps and their lengths are. -/
theorem ext' (hG : s.G.twin = t.G.twin) (hℓ : s.ℓ = t.ℓ) : s = t := by
  obtain ⟨⟨tw, htw⟩, ℓ, hm⟩ := s
  obtain ⟨⟨tw', htw'⟩, ℓ', hm'⟩ := t
  simp only at hG hℓ
  subst hG hℓ
  rfl

variable [DecidableEq F]

/-- `t` is `s` with the edge of `h₀` flipped, `h₀` having a **strictly convex diamond**: the
flips that keep a geodesic triangulation of the same surface. -/
def ConvexFlip (s t : MetricGluing F) : Prop :=
  ∃ h₀, s.G.Flippable h₀ ∧ s.G.ConvexAt s.ℓ h₀ ∧ t.G = s.G.flip h₀ ∧ t.ℓ = s.G.flipLen s.ℓ h₀

/-- **One Delaunay flip** `s ⟶ t`: `t` is `s` with an edge flipped that has two distinct faces
and is not locally Delaunay (Indermitte et al. 2001, Sec. 4; B–S Sec. 3). -/
def DelaunayFlip (s t : MetricGluing F) : Prop :=
  ∃ h₀, s.G.DelaunayFlippable s.ℓ h₀ ∧ t.G = s.G.flip h₀ ∧ t.ℓ = s.G.flipLen s.ℓ h₀

/-- The state after the Delaunay flip of `h₀` (B–S Prop. 11: it is again a metric). -/
noncomputable def flipAt (s : MetricGluing F) (h₀ : Halfedge F)
    (hd : s.G.DelaunayFlippable s.ℓ h₀) : MetricGluing F :=
  ⟨s.G.flip h₀, s.G.flipLen s.ℓ h₀, s.isMetric.flip_of_not_locallyDelaunay hd.1 hd.2⟩

theorem delaunayFlip_flipAt (s : MetricGluing F) (h₀ : Halfedge F)
    (hd : s.G.DelaunayFlippable s.ℓ h₀) : DelaunayFlip s (s.flipAt h₀ hd) :=
  ⟨h₀, hd, rfl, rfl⟩

/-- Every Delaunay flip is a flip of a strictly convex diamond (B–S Prop. 11). -/
theorem DelaunayFlip.convexFlip (h : DelaunayFlip s t) : ConvexFlip s t := by
  obtain ⟨h₀, ⟨hf, hnd⟩, hG, hℓ⟩ := h
  obtain ⟨hk, hl⟩ := s.isMetric.diamond_isTriangle h₀
  exact ⟨h₀, hf, isConvexDiamond_of_not_isLocallyDelaunay hk hl hnd, hG, hℓ⟩

theorem convexReach_of_reach (h : ReflTransGen DelaunayFlip s t) : ReflTransGen ConvexFlip s t :=
  Relation.ReflTransGen.mono (r := DelaunayFlip) (p := ConvexFlip)
    (fun _ _ h => DelaunayFlip.convexFlip h) s t h

/-- **The Delaunay states are exactly the states without a Delaunay flip** (with the self-glued
edges of `isDelaunay_iff`): the end of every maximal flip sequence is Delaunay.

[formula: it-delaunay-flip-algorithm] -/
theorem isDelaunay_iff_forall_not_delaunayFlip (s : MetricGluing F) :
    s.G.IsDelaunay s.ℓ ↔ ∀ t, ¬ DelaunayFlip s t := by
  rw [s.isMetric.isDelaunay_iff]
  constructor
  · rintro hno t ⟨h₀, hd, -, -⟩
    exact hno h₀ hd
  · intro hno h₀ hd
    exact hno _ (s.delaunayFlip_flipAt h₀ hd)

/-! ### Invariants of flip sequences -/

/-- **χ is unchanged** along any sequence of flips. -/
theorem eulerChar_eq_of_reach [Fintype F] (h : ReflTransGen ConvexFlip s t) :
    t.G.eulerChar = s.G.eulerChar := by
  induction h with
  | refl => rfl
  | tail _ hst ih =>
    obtain ⟨h₀, hf, -, hG, -⟩ := hst
    rw [hG, flip_eulerChar hf, ih]

/-- **The number of boundary halfedges is unchanged.** -/
theorem numBoundary_eq_of_reach [Fintype F] (h : ReflTransGen ConvexFlip s t) :
    t.G.numBoundary = s.G.numBoundary := by
  induction h with
  | refl => rfl
  | tail _ hst ih =>
    obtain ⟨h₀, -, -, hG, -⟩ := hst
    rw [hG, flip_numBoundary h₀, ih]

/-- **The boundary is unchanged, with its lengths**: a relabelling `e` of the halfedges maps the
boundary halfedges onto the boundary halfedges, each with its length. -/
theorem boundary_equiv_of_reach (h : ReflTransGen ConvexFlip s t) :
    ∃ e : Halfedge F ≃ Halfedge F, ∀ x,
      (t.G.twin (e x) = e x ↔ s.G.twin x = x) ∧ (s.G.twin x = x → t.ℓ (e x) = s.ℓ x) := by
  induction h with
  | refl => exact ⟨Equiv.refl _, fun x => ⟨Iff.rfl, fun _ => rfl⟩⟩
  | @tail b c _ hbc ih =>
    obtain ⟨e, he⟩ := ih
    obtain ⟨h₀, hf, -, hG, hℓ⟩ := hbc
    refine ⟨e.trans (b.G.flipPerm h₀), fun x => ⟨?_, fun hx => ?_⟩⟩
    · rw [Equiv.trans_apply, hG, flip_isBoundary_iff]
      exact (he x).1
    · have hb : b.G.twin (e x) = e x := (he x).1.2 hx
      rw [Equiv.trans_apply, hℓ, flipLen_flipPerm hf]
      split_ifs with h1
      · rcases h1 with h1 | h1
        · rw [h1] at hb
          exact absurd hb hf.twin_ne
        · rw [h1, twin_twin] at hb
          exact absurd hb.symm hf.twin_ne
      · exact (he x).2 hx

/-- **The total area is unchanged.** -/
theorem totalArea_eq_of_reach [Fintype F] (h : ReflTransGen ConvexFlip s t) :
    totalArea t.ℓ = totalArea s.ℓ := by
  induction h with
  | refl => rfl
  | @tail b c _ hbc ih =>
    obtain ⟨h₀, hf, hc, -, hℓ⟩ := hbc
    rw [hℓ, b.isMetric.totalArea_flip hf hc, ih]

/-- **Every vertex keeps its cone angle**: the vertices of `s` and `t` are in bijection, with
the same cone angles (layer 2's `coneAngle_flip` along the sequence). -/
theorem coneAngle_equiv_of_reach [Fintype F] (h : ReflTransGen ConvexFlip s t) :
    ∃ e : s.G.Vertex ≃ t.G.Vertex, ∀ v, t.G.coneAngle t.ℓ (e v) = s.G.coneAngle s.ℓ v := by
  induction h with
  | refl => exact ⟨Equiv.refl _, fun _ => rfl⟩
  | @tail b c _ hbc ih =>
    obtain ⟨e, he⟩ := ih
    obtain ⟨h₀, hf, hc, hG, hℓ⟩ := hbc
    obtain ⟨G', ℓ', hm'⟩ := c
    simp only at hG hℓ
    subst hG hℓ
    exact ⟨e.trans (flipVertexEquiv hf), fun v => by
      rw [Equiv.trans_apply]
      exact (b.isMetric.coneAngle_flip hf hc (e v)).trans (he v)⟩

/-! ### The potential and termination -/

variable [Fintype F]

/-- The potential of a state: the global lifting energy `Σ_faces A (a² + b² + c²)`. -/
noncomputable def potential (s : MetricGluing F) : ℝ := liftPotential s.ℓ

/-- **A Delaunay flip lowers the potential strictly.**

[formula: it-delaunay-potential] -/
theorem DelaunayFlip.potential_lt (h : DelaunayFlip s t) : t.potential < s.potential := by
  obtain ⟨h₀, ⟨hf, hnd⟩, -, hℓ⟩ := h
  unfold potential
  rw [hℓ]
  exact s.isMetric.liftPotential_flip_lt hf hnd

/-- A Delaunay flip lowers the harmonic index strictly (B–S Prop. 12). -/
theorem DelaunayFlip.harmonicPotential_lt (h : DelaunayFlip s t) :
    harmonicPotential t.ℓ < harmonicPotential s.ℓ := by
  obtain ⟨h₀, ⟨hf, hnd⟩, -, hℓ⟩ := h
  rw [hℓ]
  exact s.isMetric.harmonicPotential_flip_lt hf hnd

theorem potential_lt_of_transGen (h : TransGen DelaunayFlip s t) : t.potential < s.potential := by
  induction h with
  | single h => exact h.potential_lt
  | tail _ h ih => exact h.potential_lt.trans ih

/-- **No state repeats** in a sequence of Delaunay flips (no hypothesis needed): the algorithm
never cycles. -/
theorem DelaunayFlip.transGen_ne (h : TransGen DelaunayFlip s t) : s ≠ t := by
  rintro rfl
  exact lt_irrefl _ (potential_lt_of_transGen h)

theorem potential_le_of_reach (h : ReflTransGen DelaunayFlip s t) : t.potential ≤ s.potential := by
  rcases h.cases_head with rfl | ⟨u, hu, hut⟩
  · exact le_rfl
  · exact (potential_lt_of_transGen (TransGen.head' hu hut)).le

theorem harmonicPotential_le_of_reach (h : ReflTransGen DelaunayFlip s t) :
    harmonicPotential t.ℓ ≤ harmonicPotential s.ℓ := by
  induction h with
  | refl => exact le_rfl
  | tail _ h ih => exact h.harmonicPotential_lt.le.trans ih

/-- The states reachable from `s` by Delaunay flips. -/
def reach (s : MetricGluing F) : Set (MetricGluing F) := {t | ReflTransGen DelaunayFlip s t}

/-- **Termination, conditional**: if only finitely many states are reachable from `s` by
Delaunay flips, every sequence of Delaunay flips from `s` is finite (`s` is accessible for the
converse of the flip relation). The proof: the reachable states below the potential of a state
are a finite set that shrinks strictly with every flip.

[formula: it-delaunay-flip-algorithm] -/
theorem acc_of_finite (hfin : (reach s).Finite) : Acc (Function.swap DelaunayFlip) s := by
  classical
  set R := hfin.toFinset
  let rank : MetricGluing F → ℕ := fun t => (R.filter fun u => u.potential < t.potential).card
  have step : ∀ t u, t ∈ reach s → DelaunayFlip t u → u ∈ reach s ∧ rank u < rank t := by
    intro t u ht htu
    have hu : u ∈ reach s := ReflTransGen.tail ht htu
    have hlt := htu.potential_lt
    refine ⟨hu, Finset.card_lt_card ((Finset.ssubset_iff_of_subset ?_).2 ⟨u, ?_, ?_⟩)⟩
    · intro w
      simp only [Finset.mem_filter]
      exact fun ⟨hw, h⟩ => ⟨hw, h.trans hlt⟩
    · simp only [Finset.mem_filter, R, Set.Finite.mem_toFinset]
      exact ⟨hu, hlt⟩
    · simp
  suffices ∀ n, ∀ t ∈ reach s, rank t = n → Acc (Function.swap DelaunayFlip) t from
    this _ s ReflTransGen.refl rfl
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro t ht hn
    refine Acc.intro t fun u htu => ?_
    obtain ⟨hu, hlt⟩ := step t u ht htu
    exact ih _ (hn ▸ hlt) u hu rfl

/-- If every flip class is finite, the converse of the flip relation is well founded. -/
theorem wellFounded_of_finite (h : ∀ s : MetricGluing F, (reach s).Finite) :
    WellFounded (Function.swap (DelaunayFlip (F := F))) :=
  ⟨fun s => acc_of_finite (h s)⟩

omit [Fintype F] in
/-- **There is no infinite sequence of Delaunay flips** from an accessible state. -/
theorem not_exists_infinite_flips (hacc : Acc (Function.swap DelaunayFlip) s) :
    ¬ ∃ seq : ℕ → MetricGluing F, seq 0 = s ∧ ∀ n, DelaunayFlip (seq n) (seq (n + 1)) := by
  induction hacc with
  | intro t _ ih =>
    rintro ⟨seq, h0, hs⟩
    exact ih (seq 1) (by rw [← h0]; exact hs 0) ⟨fun n => seq (n + 1), rfl, fun n => hs (n + 1)⟩

omit [Fintype F] in
/-- **A Delaunay state is reachable** from an accessible state (Indermitte et al., Thm. 1; B–S
Thm. 15 / Prop. 12: existence of a Delaunay triangulation by the flip algorithm). -/
theorem exists_delaunay_of_acc (hacc : Acc (Function.swap DelaunayFlip) s) :
    ∃ t, ReflTransGen DelaunayFlip s t ∧ t.G.IsDelaunay t.ℓ := by
  induction hacc with
  | intro t _ ih =>
    by_cases hD : t.G.IsDelaunay t.ℓ
    · exact ⟨t, ReflTransGen.refl, hD⟩
    · obtain ⟨u, htu⟩ : ∃ u, DelaunayFlip t u := by
        by_contra hno
        push Not at hno
        exact hD ((isDelaunay_iff_forall_not_delaunayFlip t).2 hno)
      obtain ⟨w, hw, hD'⟩ := ih u htu
      exact ⟨w, ReflTransGen.head htu hw, hD'⟩

/-- **A bound on the number of flips**: a sequence of `n` Delaunay flips from `s` visits `n + 1`
distinct reachable states, so `n + 1 ≤ |reach s|`. -/
theorem flips_length_le (hfin : (reach s).Finite) (seq : ℕ → MetricGluing F) (n : ℕ)
    (h0 : seq 0 = s) (hs : ∀ i < n, DelaunayFlip (seq i) (seq (i + 1))) :
    n + 1 ≤ hfin.toFinset.card := by
  have mono : ∀ i k, i + k + 1 ≤ n → (seq (i + k + 1)).potential < (seq i).potential := by
    intro i k
    induction k with
    | zero => intro hk; exact (hs i (by omega)).potential_lt
    | succ k ih =>
      intro hk
      have h1 := (hs (i + k + 1) (by omega)).potential_lt
      have h2 := ih (by omega)
      rw [show i + (k + 1) + 1 = i + k + 1 + 1 by ring]
      exact h1.trans h2
  have mem : ∀ i ≤ n, seq i ∈ reach s := by
    intro i
    induction i with
    | zero => intro _; rw [h0]; exact ReflTransGen.refl
    | succ i ih => intro hi; exact ReflTransGen.tail (ih (by omega)) (hs i (by omega))
  have ne : ∀ i j, i < j → j ≤ n → seq i ≠ seq j := by
    intro i j hij hj e
    have := mono i (j - i - 1) (by omega)
    rw [show i + (j - i - 1) + 1 = j by omega, e] at this
    exact lt_irrefl _ this
  have := Finset.card_le_card_of_injOn seq (s := Finset.range (n + 1)) (t := hfin.toFinset)
    (fun i hi => by
      have hi' : i < n + 1 := by simpa using hi
      simpa using mem i (by omega))
    (by
      intro i hi j hj e
      have hi : i < n + 1 := by simpa using hi
      have hj : j < n + 1 := by simpa using hj
      by_contra hne
      rcases Nat.lt_or_gt_of_ne hne with h | h
      · exact ne i j h (by omega) e
      · exact ne j i h (by omega) e.symm)
  simpa using this

/-- **The flip algorithm, conditional on finiteness** (Indermitte et al. 2001, Thm. 1;
Bobenko–Springborn 2007, Prop. 11, 12, Thm. 15): if the states reachable from `s` by Delaunay
flips are finitely many, then
1. every sequence of Delaunay flips from `s` is finite (`Acc`), with at most `|reach s| − 1`
   flips (`flips_length_le`);
2. some Delaunay state is reachable;
3. every reachable state `t` without a further Delaunay flip (the end of a maximal sequence) is
   Delaunay, with the Euler characteristic, the number of boundary halfedges, the total area and
   the cone angles of `s`.

[formula: it-delaunay-flip-algorithm] -/
theorem flip_algorithm (hfin : (reach s).Finite) :
    Acc (Function.swap DelaunayFlip) s ∧
      (∃ t, ReflTransGen DelaunayFlip s t ∧ t.G.IsDelaunay t.ℓ) ∧
      ∀ t, ReflTransGen DelaunayFlip s t → (∀ u, ¬ DelaunayFlip t u) →
        t.G.IsDelaunay t.ℓ ∧ t.G.eulerChar = s.G.eulerChar ∧
          t.G.numBoundary = s.G.numBoundary ∧ totalArea t.ℓ = totalArea s.ℓ ∧
          ∃ e : s.G.Vertex ≃ t.G.Vertex, ∀ v, t.G.coneAngle t.ℓ (e v) = s.G.coneAngle s.ℓ v := by
  have hacc := acc_of_finite hfin
  refine ⟨hacc, exists_delaunay_of_acc hacc, fun t ht hmax => ?_⟩
  have hc := convexReach_of_reach ht
  exact ⟨(isDelaunay_iff_forall_not_delaunayFlip t).2 hmax, eulerChar_eq_of_reach hc,
    numBoundary_eq_of_reach hc, totalArea_eq_of_reach hc, coneAngle_equiv_of_reach hc⟩

/-! ## The finiteness hypothesis in the form of the literature -/

/-- **Isomorphic states**: an isomorphism of the gluing data (a relabelling of the faces with a
rotation of each, `Gluing.Iso`) that carries the lengths. Isomorphic states are the same
triangulated surface. -/
def IsIso (s t : MetricGluing F) : Prop :=
  ∃ ψ : Gluing.Iso s.G t.G, ∀ h, t.ℓ (ψ.toHalfedge h) = s.ℓ h

/-- The relabelling of the halfedges given by a face bijection and a rotation of each face. -/
def relabel (p : (F ≃ F) × (F → Fin 3)) : Halfedge F ≃ Halfedge F where
  toFun h := (p.1 h.1, h.2 + p.2 h.1)
  invFun h := (p.1.symm h.1, h.2 - p.2 (p.1.symm h.1))
  left_inv h := by simp
  right_inv h := by simp

omit [DecidableEq F] in
/-- **A state has finitely many isomorphic states** (with the same faces `F`): an isomorphism is
determined by a face bijection and a rotation of each face. -/
theorem finite_isIso (s : MetricGluing F) : {t | IsIso s t}.Finite := by
  classical
  let g : (F ≃ F) × (F → Fin 3) → (Halfedge F → Halfedge F) × (Halfedge F → ℝ) := fun p =>
    (fun y => relabel p (s.G.twin ((relabel p).symm y)), fun y => s.ℓ ((relabel p).symm y))
  let κ : MetricGluing F → (Halfedge F → Halfedge F) × (Halfedge F → ℝ) := fun t => (t.G.twin, t.ℓ)
  have hκ : Function.Injective κ := fun t u e =>
    ext' (congrArg Prod.fst e) (congrArg Prod.snd e)
  refine ((Set.finite_range g).preimage hκ.injOn).subset ?_
  rintro t ⟨ψ, hψ⟩
  refine ⟨(ψ.faces, ψ.shift), ?_⟩
  have hrel : ∀ y, relabel (ψ.faces, ψ.shift) y = ψ.toHalfedge y := fun _ => rfl
  have hrel' : ∀ y, (relabel (ψ.faces, ψ.shift)).symm y = ψ.toHalfedge.symm y := fun _ => rfl
  simp only [g, κ, Prod.mk.injEq]
  constructor
  · funext y
    rw [hrel', hrel, ψ.toHalfedge_twin, Equiv.apply_symm_apply]
  · funext y
    rw [hrel', ← hψ, Equiv.apply_symm_apply]

/-- The **flip class** of `s`: the states reachable from `s` by flips of strictly convex diamonds.
Each is a geodesic triangulation of the same piecewise flat surface with the same vertices
(the surface itself is layer 4, not modelled here). -/
def flipClass (s : MetricGluing F) : Set (MetricGluing F) := {t | ReflTransGen ConvexFlip s t}

/-- **The finiteness hypothesis** of the termination theorem, as the literature uses it: for every
length `L`, the triangulations in the flip class of `s` whose edges are all at most `L` long are
finitely many up to isomorphism. On a compact piecewise flat surface this holds because there
are only finitely many geodesic arcs of bounded length between the (finitely many) vertices
(B–S, proof of Thm. 15; Indermitte et al., proof of Thm. 1, in the form "takes values in a
discrete set"); it needs the metric surface and is cited (`CITED.md`,
`it-delaunay-flip-termination`). -/
def BoundedLengthFinite (s : MetricGluing F) : Prop :=
  ∀ L : ℝ, ∃ S : Finset (MetricGluing F),
    ∀ t ∈ flipClass s, (∀ h, t.ℓ h ≤ L) → ∃ u ∈ S, IsIso u t

/-- Along Delaunay flips every edge stays shorter than `√(H(s)·A(s))`: the harmonic index does not
increase and the area is constant (B–S's argument for bounded edge lengths). -/
theorem length_le_of_reach (h : ReflTransGen DelaunayFlip s t) (x : Halfedge F) :
    t.ℓ x ≤ √(harmonicPotential s.ℓ * totalArea s.ℓ) := by
  have h1 := t.isMetric.sq_le_harmonicPotential_mul_totalArea x
  have h2 := harmonicPotential_le_of_reach h
  have h3 := totalArea_eq_of_reach (convexReach_of_reach h)
  have hA : 0 ≤ totalArea t.ℓ := Finset.sum_nonneg fun f _ => (triArea_pos (t.isMetric.tri f)).le
  have : t.ℓ x ^ 2 ≤ harmonicPotential s.ℓ * totalArea s.ℓ := by
    rw [← h3]
    exact h1.trans (mul_le_mul_of_nonneg_right h2 hA)
  exact (le_abs_self _).trans (Real.abs_le_sqrt this)

/-- **The literature's finiteness hypothesis gives a finite reachable set**, hence termination
(`flip_algorithm`).

[formula: it-delaunay-flip-algorithm] -/
theorem reach_finite_of_boundedLengthFinite (hB : BoundedLengthFinite s) : (reach s).Finite := by
  obtain ⟨S, hS⟩ := hB (√(harmonicPotential s.ℓ * totalArea s.ℓ))
  refine (S.finite_toSet.biUnion fun u _ => finite_isIso u).subset fun t ht => ?_
  obtain ⟨u, hu, hut⟩ := hS t (convexReach_of_reach ht)
    (length_le_of_reach ht)
  exact Set.mem_biUnion hu hut

end MetricGluing

end DDG.Delaunay
