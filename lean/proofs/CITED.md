# Cited, not machine-checked

Results the library relies on, or that the paper uses to motivate it, which are **not** proved in
`lean/proofs`. Each entry has an id (the same `[formula: id]` tag the Rust code may carry;
`scripts/check-formulas.sh` reads the ids from the headings below), an exact reference, the reason
it is not formalized, and whether the correctness of the library *relative to the canonical
discrete problem of `CLAUDE.md`* depends on it.

"Paper" is Soliman, Slepčev, Crane, *Optimal Cone Singularities for Conformal Flattening*,
ACM TOG 37(4), 2018. "Thesis" is Soliman, *Optimal Cone Singularities for Conformal Flattening*,
MS thesis, CMU 2018. "BFF" is Sawhney, Crane, *Boundary First Flattening*, ACM TOG 37(1), 2017
(page and equation numbers of the authors' version,
`https://www.cs.cmu.edu/~kmcrane/Projects/BoundaryFirstFlattening/paper.pdf`). Mathlib is the
revision pinned in `lake-manifest.json` (tag `v4.34.1`); a search of that revision finds no
subdifferential, no convex conjugate / Fenchel duality, no Wasserstein distance, no
Laplace–Beltrami operator or Gaussian curvature of a Riemannian surface, no classification of
surfaces, and Sobolev spaces only in the form of the Euclidean Sobolev inequality and Bessel
potential spaces.

## What the machine-checked part covers, and what it does not

Every theorem of `lean/proofs` is a statement in **exact real arithmetic** about the formulas the
library evaluates. Nothing is proved about floating-point rounding, and the Rust code is not
connected to the Lean statements mechanically (the correspondence is by reading, by the
`[formula: id]` tags, and by differential tests).

For the cone solver, the claim "the returned `(u, µ)` solves the canonical discrete problem" has
this status:

* **If the optimality system holds exactly** (`L u = µ − Ω`, `Lᵀφ = −wE m u`, the inclusions),
  `(u, µ)` is a global minimizer: machine-checked (`optimality-system-sufficient`). Likewise the
  duality-gap certificate for an exactly feasible primal point and an exactly feasible dual point
  (`duality_gap_certificate`, `adjoint_gap_certificate`).
* **A computed triple satisfies none of these equations exactly.** What is machine-checked for
  arbitrary `(u, φ, µ)` (`DdgProofs/Certificate.lean`, ids `duality-gap-certificate`,
  `duality-gap-inexact`), with `r = L u − µ + Ω`, `φ̂` the clipped `φ`,
  `D = dualObjectiveOfAdjoint φ̂`:
  - `D ≤` optimum, provided `(Lᵀφ̂)_i = 0` wherever `wE_i m_i = 0` (no condition if all
    `wE_i m_i > 0`);
  - if `µ_i + r_i` is admissible at every vertex: `(u, µ + r)` is exactly feasible,
    optimum `≤ J(u, µ + r) ≤ J(u, µ) + λ Σ wR_i |r_i|`, and `(u, µ + r)` is suboptimal by at most
    `J(u, µ) − D + λ Σ wR_i |r_i|`; every minimizer `u*` has
    `¼ Σ wE m (u − u*)² ≤` that width;
  - if `µ_i + r_i` is **not** admissible somewhere (a forbidden or prescribed vertex with
    `r_i ≠ 0`; a vertex at a bound whose residual points outward), nothing is proved about the
    original problem. Proved instead: `(u, µ)` is within the complementarity gap of the optimum
    of the problem with angle defects `Ω − r` (`backward_certificate`), and the repaired point is
    certified for the problem with defects changed only at the offending vertices
    (`repaired_certificate`). The examples in `Certificate.lean` show that the forward statement
    is false without admissibility.
  - `J(u, µ) − D` alone (the computed "duality gap") bounds nothing: it can be negative.
* **Not covered at all:** rounding in evaluating `L u`, `r`, `D`, `J` and the tests
  (`µ_i + r_i` admissible, `φ̂_i` in the domain); that the Newton iteration terminates
  (`semismooth-newton-convergence`, cited); existence of a minimizer
  (`discrete-minimizer-existence`, cited).

The cited results below are not needed for any of the machine-checked statements. They explain
why the canonical discrete problem is the right problem (smooth theory), why the solver can be
expected to terminate (Newton convergence), what the optional post-processing does (cone merging,
genus reporting), and — for the flattener — the smooth theory behind Boundary First Flattening
and three topological facts (classification of surfaces, tree–cotree cuts, injectivity from the
boundary).

---

## Smooth formulation

### [formula: smooth-yamabe-equation] Yamabe equation with cone singularities

- **Statement.** Under `g̃ = e^{2u} g` the Gaussian curvature changes by `Δu = K − e^{2u} K̃`, with
  the Cherrier boundary condition `∂u/∂n = κ − e^{u} κ̃`; for a polyhedral cone metric
  `Δu = K − Σ θ_i δ_{p_i}`.
- **Reference.** Paper Sec. 3.1, Eq. 1; Aubin, *Some Nonlinear Problems in Riemannian Geometry*
  (1998); Troyanov 1989, Def. 1 (cone metrics).
- **Why not formalized.** Mathlib has no Gaussian curvature, Laplace–Beltrami operator or
  conformal-change formulas for Riemannian surfaces.
- **Does correctness depend on it?** No. It is the modelling step that leads to the constraint
  `L u = µ − Ω`. Its discrete counterpart **is** proved:
  `TriMesh.hasDerivAt_angleDefect_conformalScale` (`[formula: discrete-yamabe-linearization]`).

### [formula: cone-problem-relaxation] Relaxation to measures and existence of minimizers

- **Statement.** The cone placement problem (Problem 5) is relaxed to a convex problem over
  Radon measures (Problem 9), which has a minimizer; the `L¹` version (Problem 8) in general
  has none.
- **Reference.** Paper Sec. 4.1–4.2, Problems 5, 7, 8, 9; Thesis Ch. 3.4 and Ch. 4.3.
- **Why not formalized.** Needs `M(M) = C(M)*` (Riesz–Markov), weak-* compactness and the direct
  method on a Riemannian surface.
- **Does correctness depend on it?** No. The canonical discrete problem is taken as the ground
  truth.

### [formula: smooth-predual-duality] Pre-dual problem and smooth optimality system

- **Statement.** Fenchel–Rockafellar duality between Problem 10 (pre-dual, functions) and
  Problem 9 (measures); the optimality system `Δu = K − µ`, `Δφ = u`, `µ ∈ ∂1_λ(φ)` (System 11)
  and its generalization with weights, regions and prescribed sites.
- **Reference.** Paper Sec. 3.4 (Eqs. 2–4), Sec. 4.3 (Problem 10), Sec. 4.4 (System 11),
  App. A.2, App. A.3 (Prop. A.3, Problems A.1/A.2, generalized optimality conditions); Brezis,
  *Functional Analysis, Sobolev Spaces and PDE* (2010); Thesis Ch. 4.1, 4.3.
- **Why not formalized.** Mathlib has no convex conjugate, subdifferential or Fenchel–Rockafellar
  theorem, and the duality is between `C(M)` and `M(M)`.
- **Does correctness depend on it?** No. The finite-dimensional analogues are proved directly:
  sufficiency and necessity of the discrete optimality system
  (`[formula: optimality-system-sufficient]`, `[formula: optimality-system-necessary]`), weak
  duality, the gap identity and strong duality at a solution
  (`[formula: duality-gap-certificate]`).

### [formula: poisson-measure-data] Poisson equation with measure data

- **Statement.** Existence of the Green's function with a logarithmic bound; for a signed Radon
  measure `µ` the Dirichlet problem `Δu = µ` has a unique solution in `W₀^{1,p}`, `1 < p < 2`, with
  `‖∇u‖_{L^p} ≤ c‖µ‖_M`; hence a solution operator `S : M(M) → L²(M)` which is the adjoint of the
  solution operator `S* : L²(M) → C(M)`; extension to boundary data in `H^{3/2}(∂M)`.
- **Reference.** Paper App. A.1, Thm. A.1 (Aubin 1998, Thm. 4.17), Thm. A.2 (Littman,
  Stampacchia, Weinberger 1963, Thm. 5.1), App. A.1.1; Thesis Ch. 4.2.
- **Why not formalized.** Elliptic regularity and Sobolev embeddings on manifolds are not in
  Mathlib.
- **Does correctness depend on it?** No. Discretely the state equation is a linear system whose
  matrix is positive definite (`[formula: cotan-laplacian-interior-posdef]`, proved) or, on a
  closed mesh, positive semidefinite with kernel the constants
  (`[formula: cotan-laplacian-kernel]`, proved).

### [formula: optimal-measure-support] Support of the optimal measure

- **Statement.** `supp µ⁺ ⊆ {φ = +λ w_R}`, `supp µ⁻ ⊆ {φ = −λ w_R}`.
- **Reference.** Paper Prop. 4.1, Prop. C.1.
- **Why not formalized.** Smooth statement about Radon measures (Jordan decomposition,
  Radon–Nikodym).
- **Does correctness depend on it?** No. The discrete statement is proved
  (`eq_softThreshold_add_iff_cases`, `[formula: soft-threshold-fixed-point]`): `φ_i = λ w_R` where
  `µ_i > 0`, `φ_i = −λ w_R` where `µ_i < 0`, `|φ_i| ≤ λ w_R` everywhere.

### [formula: optimal-measure-regularity] The smooth minimizer has no isolated atoms

- **Statement.** `‖u‖_{L^∞}` bound for the optimal scale factor; the optimizing measure lies in
  `H₀¹(M)*`, so it contains no Dirac masses.
- **Reference.** Paper App. C, Lemma C.2, Lemma C.3, Thm. C.4 (using Saff, Totik 1997,
  Cor. 3.3); Thesis Ch. 4.4.
- **Why not formalized.** Potential theory and Sobolev duality.
- **Does correctness depend on it?** No. It explains why, on very fine meshes, the discrete
  minimizer may contain small clusters of cones (paper Fig. 9).

### [formula: cone-rounding-stability] Stability of rounding a measure to cones

- **Statement.** For probability measures on a ball, `‖µ − ν‖_{H^{-2}} ≲ W₂(µ, ν)`; consequently
  replacing a small cluster of cones by one cone of the same total angle changes the distortion
  by a small amount; linear combinations of Dirac deltas are weak-* dense.
- **Reference.** Paper App. C, Thm. C.5, Cor. C.6, Lemma C.7 (Mainini 2012, Prop. 3.8; Helms 2009,
  Thm. 3.2.12); used in paper Sec. 6.1 ("Extracting cones").
- **Why not formalized.** Mathlib has no optimal transport / Wasserstein distance and no Green's
  function estimates.
- **Does correctness depend on it?** Only for the optional post-processing of paper Sec. 6.1
  (merging edge-connected clusters of cones into one cone). A merged configuration is feasible but
  is no longer the minimizer; that it is nearly as good is exactly this cited estimate. The
  optimality certificate refers to the configuration *before* merging and does not depend on it.

### [formula: moreau-yosida-convergence] Convergence of the regularized problems

- **Statement.** The minimizers `µ_γ` of the regularized problems (Problem 12, System 13)
  converge weak-* to the minimizer as `γ → 0`.
- **Reference.** Paper Sec. 4.6, Eqs. 12–13; Thesis Ch. 4.5 (Sec. 4.5.2).
- **Why not formalized.** Γ-convergence / weak-* compactness argument in the smooth setting; the
  finite-dimensional analogue was not attempted.
- **Does correctness depend on it?** No. Stage I only produces a starting point for Stage II,
  whose result is certified independently. Proved discretely: a zero of the regularized residual
  minimizes the regularized objective (`[formula: regularized-optimality-sufficient]`), and its
  objective exceeds the optimum by at most `Σ_i (γ_i/2) µ*_i²`
  (`[formula: regularized-suboptimality]`).

---

## Discretization

### [formula: variational-discretization] Variational discretization of the measure space

- **Statement.** With the state discretized by piecewise linear finite elements and the measure
  left undiscretized, the semidiscrete problem (Problem B.1) has a minimizer that is a linear
  combination of Dirac measures at the vertices, unique among such combinations for a given
  state; the fully discrete optimality systems are Systems B.2 and B.3.
- **Reference.** Paper App. B.1; Casas, Clason, Kunisch 2012, Thm. 3.2; Hinze 2005; Thesis Ch. 5.3.
- **Why not formalized.** Finite element theory over measure spaces.
- **Does correctness depend on it?** No. This result is what relates the canonical discrete
  problem to the smooth one; the canonical problem itself is fully discrete and is treated
  directly.

### [formula: cotan-fem-convergence] The cotan finite element method converges

- **Statement.** `L u = M f` discretizes `Δu = f` on a surface, and its solutions converge to
  those of the Laplace–Beltrami problem under refinement.
- **Reference.** Paper Sec. 5.1.3; MacNeal 1949; Pinkall, Polthier 1993.
- **Why not formalized.** Convergence of finite elements on surfaces needs Sobolev spaces on a
  surface, which Mathlib does not have. Proved instead (`DdgProofs/TriangleIntegral/Fem.lean`,
  `[formula: cotan-fem-stiffness]`): `L` is the stiffness matrix of piecewise linear elements,
  as Lebesgue integrals over isometric layouts of the faces — `∫_T ∇λᵢ·∇λⱼ = −½ cot θ_k`,
  `∫_T |∇f|² = ½ Σ cot θ_k (f_i − f_j)²`, `Σ_f ∫_{T_f} |∇u|² = uᵀLu` — and the lumped mass is the
  row sum of the Galerkin mass `∫_T λᵢ λⱼ` (`[formula: app-fem-galerkin-mass]`).
- **Does correctness depend on it?** No.

### [formula: semismooth-newton-convergence] Convergence of the semismooth Newton method

- **Statement.** `P_λ` is semismooth with generalized derivative `A⁺ + A⁻`; the semismooth Newton
  iteration converges locally superlinearly, with mesh-independent rate, and the active-set
  stopping criterion is justified.
- **Reference.** Paper Sec. 5.3 (Eqs. 23–26), Sec. 5.3.3; Hintermüller, Ito, Kunisch 2002,
  Lemma 3.1; Ulbrich 2011, Ch. 5 and 7; Clason, Schiela 2017; Thesis Ch. 5.2 (Prop. 60).
- **Why not formalized.** Mathlib has no Newton (semismooth) differentiability; the convergence
  analysis was not attempted.
- **Does correctness depend on it?** No. What is proved is what the stopping rule needs
  (`[formula: newton-stopping-rule]`, `[formula: soft-threshold-active-sets]`): the residual is
  piecewise affine, a step that keeps the active sets lands on an exact zero, and an exact zero
  is a global minimizer. That the iteration reaches such a step in finitely many iterations is
  **not** proved; a run that does not is detected at run time.

### [formula: discrete-minimizer-existence] Existence of a minimizer of the canonical problem

- **Statement.** The canonical discrete problem has a minimizer (a convex quadratic plus a
  polyhedral function that is bounded below on a polyhedron attains its minimum).
- **Reference.** Frank, Wolfe 1956 (quadratic programs); Rockafellar, *Convex Analysis* (1970),
  Cor. 27.3.1; elementary (Weierstrass) when all `wE_i m_i > 0`.
- **Why not formalized.** Not attempted; not needed for any statement about solver output.
- **Does correctness depend on it?** No. Every theorem here is conditional on a given candidate
  `(u, φ, µ)`; uniqueness is proved, existence is not.

### [formula: optimality-system-necessary-closed] Necessity of the optimality system on closed meshes

- **Statement.** On a closed mesh (`L` positive semidefinite with kernel the constants) every
  minimizer of the canonical problem admits an adjoint state `φ` solving the optimality system
  (the additive constant of `φ` is the multiplier of the constraint `Σµ = ΣΩ`).
- **Reference.** Rockafellar, *Convex Analysis* (1970), Thm. 28.2–28.3 and Cor. 28.2.2
  (polyhedral case; no constraint qualification needed).
- **Why not formalized.** Proved here only for invertible `L` (Dirichlet case,
  `IsMinimizer.existsUnique_kkt`). The closed case needs, in addition, the range of the singular
  symmetric matrix and a common-point argument for the per-vertex subdifferential intervals; it
  was not done.
- **Does correctness depend on it?** No. Soundness of the solver output uses sufficiency, which
  is proved for every `L`. Necessity is only the converse ("every minimizer is found by solving
  the optimality system").
- **Also covers** the closed-mesh converse of `reduced-optimality-conditions` (factor-once method):
  that every minimizer of the reduced problem with the sum constraint admits a multiplier `ν` with
  `−(y_j + ν) ∈ ∂g_j(µ_j)`. The Dirichlet converse is proved
  (`reducedOptimal_of_isMinimizer_dirichlet`); the direction the solver relies on (reduced conditions
  ⟹ `KKT` ⟹ minimizer) is proved in both cases.

---

## Topology

### [formula: genus-from-euler] Genus from the Euler characteristic

- **Statement.** A compact connected orientable surface with `b` boundary components has
  `χ = 2 − 2g − b`; summed over `C` components, `g = (2C − b − χ)/2`.
- **Reference.** Classification of compact surfaces, e.g. Massey, *Algebraic Topology: An
  Introduction* (1967), Ch. 1; used by `TriMesh::genus` in `crates/ddg-core`.
- **Why not formalized.** Mathlib has no classification of surfaces; here this is in effect the
  definition of the reported genus.
- **Does correctness depend on it?** No. The genus is reported for information; the optimization
  uses `χ` only through discrete Gauss–Bonnet, which is proved (`[formula: gauss-bonnet]`).

### [formula: tree-cotree-disk] Cutting along the complement of a dual spanning tree gives a disk

- **Statement.** (i) On a connected triangulated surface, let `T` be a set of interior edges
  that contains no cycle in the graph obtained by contracting every boundary loop to a vertex.
  Then the dual graph with the edges crossing `T` removed is connected, so it has a spanning tree
  `T*` that crosses no edge of `T`. (ii) Cutting all interior edges not crossed by `T*` glues the
  faces along a tree, and the result is a topological disk.
- **Reference.** Eppstein, *Dynamic generators of topologically embedded graphs*, SODA 2003,
  Lemma 3.1 and Sec. 3 (tree–cotree decomposition); Erickson, Har-Peled, *Optimally cutting a
  surface into a disk*, Discrete Comput. Geom. 31 (2004), Sec. 2–3; for (i) also Mohar, Thomassen,
  *Graphs on Surfaces* (2001), Sec. 4.1 (cuts of the dual are boundaries of face sets); used in
  BFF Sec. 6.4 ("cut `M₀` into a disk `M` via a cut passing through all cones", citing
  Springborn, Schröder, Pinkall 2008, Sec. 3.3). Used by `plan_pieces` in `crates/ddg-bff`.
- **Why not formalized.** Needs the topology of surfaces glued from triangles (disks under
  gluing along boundary arcs); Mathlib has no such theory. What *is* proved about cuts is the
  bookkeeping: `[formula: cut-euler-characteristic]`, `[formula: bff-cone-angle-sum]`.
- **Does correctness depend on it?** No. The planner relies on it to choose its cut, but every
  chart of the cut mesh is then checked to have Euler characteristic 1 and one boundary loop
  (`ComponentTopology::is_disk`; a failure is the error `ChartNotDisk`). That check in turn rests
  on `[formula: genus-from-euler]` (cited).

### [formula: cut-leaf-pruning] Un-cutting a dangling cut edge keeps a disk a disk

- **Statement.** If a chart is a disk and a cut edge ends at an interior vertex that has no
  other cut edge and is not a cone, gluing the two sides of that edge back together gives a disk
  again (unless the two sides are the whole boundary). Conversely, extending a cut by an edge
  that ends at such a vertex keeps a disk a disk.
- **Reference.** Erickson, Har-Peled 2004 (as above), Sec. 3 ("pruning" a cut graph: removing
  degree-one vertices that are not required); gluing two consecutive boundary edges of a polygon
  is the cancellation step `a a⁻¹` of the classification of surfaces, Massey, *Algebraic
  Topology: An Introduction* (1967), Ch. 1, Sec. 7. Used by `plan_pieces` and
  `extend_isolated_slits` in `crates/ddg-bff`.
- **Why not formalized.** As for `tree-cotree-disk`.
- **Does correctness depend on it?** No, for the same reason: the charts are checked after
  cutting. The Euler characteristic is consistent with the claim by a proved identity: such an
  un-cut removes one vertex copy and one edge copy (`[formula: cut-euler-characteristic]`).

---

## Boundary First Flattening

Proved for the flattener (see the README): existence and uniqueness of the discrete Poisson
solves, the Dirichlet-to-Neumann relation and its compatibility, solvability of the Neumann
problem, the exterior-angle sums (BFF Propositions 1 and 2), the first-order Cherrier formula,
the closing of the boundary polygon and the optimal length adjustment, the Hilbert-transform
right-hand side, the Euler characteristic after cutting, the definitions of the distortion
measures, and for the map to the sphere the stereographic projection (on the sphere, inverse,
conformal, orientation preserving), the inversions of the sphere and the derivative behind the
centering Jacobian. Cited:

### [formula: bff-target-boundary-length] Target lengths of the boundary edges

- **Statement.** `ℓ*_ij = e^{(u_i + u_j)/2} ℓ_ij`: the lengths of the discretely conformally
  equivalent metric with log scale factors `u`.
- **Reference.** BFF Sec. 5, Step II, and Alg. 1 line 9; the notion of discrete conformal
  equivalence is Luo, *Combinatorial Yamabe flow on surfaces*, Commun. Contemp. Math. 6 (2004),
  and Springborn, Schröder, Pinkall, *Conformal equivalence of triangle meshes*, ACM TOG 27(3),
  2008, Def. 2.1.
- **Why not formalized.** It is a definition, the same scaling as `TriMesh.conformalScale` in
  `DdgProofs/Yamabe.lean` (for which the first-order change of the curvatures is proved,
  `[formula: discrete-yamabe-linearization]`, `[formula: bff-cherrier-boundary-curvature]`).
  Proved about the transcription (`DdgProofs/Bff/Curve.lean`): it is positive for positive `ℓ`
  and the identity for `u = 0`.
- **Does correctness depend on it?** It defines what the flattener aims at. The layout is *not*
  claimed to have these lengths: the boundary lengths are adjusted to close the curve
  (`[formula: bff-closing-lengths]`, the change is reported), and interior lengths are whatever
  the harmonic extension gives (the distortion is measured and reported).

### [formula: bff-smooth-cherrier-formula] Conformal change of curvature with boundary

- **Statement.** For a conformal map with scale factor `e^u`: `Δu = K − e^{2u} K̃` in the
  interior and `∂u/∂n = κ − e^u κ̃` on the boundary; integrated over dual cells and dual boundary
  edges these become the linear relations `A u = Ω − Ω̃` and `h = k − k̃` between curvature
  densities.
- **Reference.** BFF Sec. 4.2, Eqs. 8–12; Cherrier, *Problèmes de Neumann non linéaires sur les
  variétés riemanniennes*, J. Funct. Anal. 57 (1984); Aubin 1998, Ch. 5.
- **Sign.** BFF Eq. 11 and the proofs of its Propositions 1 and 2 print `A u = Ω − Ω̃` with the
  positive semidefinite cotan matrix `A` of its Sec. 4.1; its Alg. 1 (lines 6 and 8) and the
  reference implementation use the source term `−Ω`, i.e. `A u = Ω̃ − Ω`. The latter is the
  convention of `CLAUDE.md` (`L u = µ − Ω`, `u` the true log scale factor) and the one proved:
  `TriMesh.hasDerivAt_angleDefect_conformalScale`. The relation `k̃ = k − h = k + (A u)_B` is the
  same in both.
- **Why not formalized.** Curvature of Riemannian surfaces with boundary is not in Mathlib.
- **Does correctness depend on it?** No. The discrete statements the code uses are proved:
  the first-order change of the vertex curvatures (`[formula: bff-cherrier-boundary-curvature]`),
  and the exact sums of the resulting angles (`[formula: bff-exterior-angle-sum]`,
  `[formula: bff-cone-angle-sum]`). That the *linearized* angles are good targets is the
  modelling assumption of BFF (its Sec. 4.2, "one does not obtain exact target curvatures if `u`
  values are used to rescale edge lengths").

### [formula: bff-disk-uniformization] Flattening to a round disk

- **Statement.** A conformal map of a disk-like surface onto the unit disk exists
  (uniformization); it has constant boundary curvature, so its curvature density is the new
  length density, and the fixed-point iteration "`k̃ ∝` current dual boundary lengths, recompute
  `u`" approaches it.
- **Reference.** BFF Sec. 6.5 (the iteration; no convergence proof is given there).
- **Why not formalized.** Uniformization is not in Mathlib; the iteration is a heuristic.
- **Does correctness depend on it?** No. Proved about the step the code takes: the prescribed
  angles are positive and sum to `2π` (`[formula: bff-uniformization-angles]`), so they are valid
  input for the flattening with prescribed exterior angles. How round the result is, is not
  claimed.

### [formula: disk-injectivity-criterion] A flattening without flipped triangles and with a simple boundary is injective

- **Statement.** A piecewise linear map of a triangulated disk to the plane in which every
  triangle keeps its orientation (no flipped or degenerate triangle), the angles around every
  boundary vertex sum to less than `2π`, and the boundary polygon does not cross itself, is
  injective.
- **Reference.** Lipman, *Bijective mappings of meshes with boundary and the degree in mesh
  processing*, SIAM J. Imaging Sci. 7(2), 2014, Thm. 1 and Sec. 3 (degree argument: an
  orientation preserving simplicial map covers each point as often as the boundary winds around
  it); the smooth predecessor is Meisters, Olech, *Locally one-to-one mappings and a classical
  theorem on schlicht functions*, Duke Math. J. 30 (1963). BFF Sec. 7.5 ("Injectivity") cites
  Floater, *One-to-one piecewise linear mappings over triangulations*, Math. Comp. 72 (2003),
  Thm. 4.1, for the convex case.
- **Why not formalized.** Needs degree theory / winding numbers of polygons in the plane, which
  Mathlib has only in pieces.
- **Does correctness depend on it?** Only the meaning of the reported flag
  `ChartReport::is_injective`. The flag is computed from counts that are defined by proved
  formulas (`[formula: signed-area-2d]`, `[formula: signed-angle-2d]`,
  `[formula: segment-proper-crossing]`); that these counts being zero implies injectivity is this
  cited theorem. As the code's documentation says, boundary pieces that coincide (the two sides
  of a cut through a flat region) are contact, not crossing, and are not detected by it; the
  sampled coverage is an independent diagnostic.

### [formula: mobius-centering-existence] A centering inversion exists and is unique up to rotation

- **Statement.** For a conformal map `f` of a genus-0 surface `M` to the unit sphere and the area
  measure of `M` carried to the sphere, there is an inversion `η_c(x) = (1 − |c|²)(x + c)/|x + c|²
  + c`, `|c| < 1`, after which the center of mass `∫ η_c(x) dA` is at the origin; the composite
  map is unique up to a rotation. The Gauss–Newton iteration `c = −J⁻¹ µ`,
  `J = 2 Σ A_f (I − C_f C_fᵀ)`, is the algorithm that finds it.
- **Reference.** A. Baden, K. Crane, M. Kazhdan, *Möbius Registration*, Computer Graphics Forum
  37(5) (SGP 2018), Sec. 4 (Eq. 1–4), Sec. 4.1 (existence and uniqueness), Sec. 6.1 and
  Algorithm 1 (the iteration, with areas `A` "from the original surface"); BFF's reference code
  cites the same algorithm in `BFF::centerMobius`.
- **Why not formalized.** Existence and uniqueness are a statement about all inversions of the
  ball (a degree / convexity argument in hyperbolic space); convergence of the iteration is not
  proved in the paper either.
- **Does correctness depend on it?** Only the claim that the result is canonical up to a
  rotation. What the code relies on is proved: each step maps the sphere to itself
  (`[formula: mobius-sphere-inversion]`), the Jacobian is the derivative of the inversion at
  `c = 0` (`[formula: mobius-centering-jacobian]`). Whether `|µ|` reached the tolerance is
  measured and reported (`CenteringReport`); a step that does not decrease `|µ|` is halved, so
  `|µ|` never increases.

### [formula: spherical-triangle-area] Signed area of a spherical triangle

- **Statement.** For unit vectors `a, b, c` the geodesic triangle they span has area `E` (its
  spherical excess) with `tan(E/2) = det(a, b, c) / (1 + a·b + b·c + c·a)`; with `atan2`, `E` is
  signed by the orientation of `a, b, c` seen from outside.
- **Reference.** A. van Oosterom, J. Strackee, *The solid angle of a plane triangle*, IEEE
  Trans. Biomed. Eng. BME-30(2), 1983, 125–126, Eq. (8).
- **Why not formalized.** Mathlib has no spherical trigonometry (areas of geodesic triangles).
- **Does correctness depend on it?** Only the reported `covering_degree` and the flag
  `SphereDistortion::is_injective`; the distortion measures and the count of flipped triangles
  (the sign of `det(a, b, c)`) do not use it.

### [formula: sphere-injectivity-criterion] A map of the sphere without flipped triangles and of degree one is injective

- **Statement.** A continuous map of a closed oriented surface of genus 0 to the unit sphere that
  maps each triangle onto a geodesic triangle with positive orientation is an orientation
  preserving branched covering; its degree is the total signed area divided by `4π`; a degree of
  1 means it is a homeomorphism.
- **Reference.** The degree of a map between closed oriented manifolds and its computation by
  counting preimages with orientation: J. Milnor, *Topology from the Differentiable Viewpoint*,
  1965, §5; branched coverings of surfaces and the Riemann–Hurwitz formula, e.g. O. Forster,
  *Lectures on Riemann Surfaces*, 1981, Thm. 4.24 (a proper holomorphic, here: orientation
  preserving locally injective away from finitely many points, map of degree 1 is bijective).
  For triangle meshes: Y. Lipman, *Bijective mappings of meshes with boundary and the degree in
  mesh processing*, SIAM J. Imaging Sci. 7(2), 2014, Sec. 3.
- **Why not formalized.** Degree theory of maps between surfaces is not in Mathlib.
- **Does correctness depend on it?** Only the meaning of the reported flag
  `SphereDistortion::is_injective`; it is computed from counts and a sum that are defined by
  `[formula: spherical-triangle-area]` and the sign of a determinant.

---

## Weights

### [formula: feature-size-weights] Curvature-adapted weights

- **Statement.** `r_i` is built from the curvature density `|Ω_i|/m_i`; `wE_i := r_i`,
  `wR_i := 1/r_i`.
- **Reference.** Paper Sec. 8.1 ("Nonuniform importance"), which cites Sharp, Crane 2018 for the
  intrinsic local feature size.
- **Discrepancy in the paper.** Sec. 8.1 prints `r_i := 1/(|Ω_i|/M_ii + ε)` and then says "this
  value will be small in flat regions and large in highly curved regions": the printed formula is
  large in flat regions. The library offers both readings and names them
  (`FeatureSizeConvention::AsDescribed`, `AsPrinted`); it applies neither on its own.
- **Why not formalized.** It is a heuristic choice of input data, not a derived formula.
- **Does correctness depend on it?** No. Any weights with `wE ≥ 0`, `wR ≥ 0` give a well-formed
  canonical problem, for which all theorems hold; for `ε > 0` and positive masses both readings
  give `r_i > 0`.

### [formula: cone-feedback-rule] The reweighting rule of the feedback loop

- **Statement.** `crates/ddg-optimal-cones/src/feedback.rs` (`Feedback::observe_probe`,
  `Feedback::observe`): first on the surface flattened without cones (the cone-free probe, where
  every component has a boundary), then after every solve and its flattening, every interior
  vertex with `wE_i > 0` that was not reweighted before and has `|δ_i| > threshold` (default 20;
  `δ` the scale mismatch, `cone-feedback-mismatch`) gets `wE_i ← reduced_weight · wE_i` (default
  `0`), and the problem is solved again, until a round after a solve picks no vertex, `max_rounds`
  solves were done, more than `max_reweighted_fraction` of the vertices would be reweighted, or a
  connected component would be left without a positive weight (then nothing of that round is
  applied; the probe never ends the loop).
- **Reference.** None: this is a **design choice, not a formula**. The paper (Sec. 7.2, Fig. 20)
  says that needle-like outliers "are ignored" because they carry almost no area, without saying
  how; the crate README ("Needle-like spikes") measured that zero distortion weight at the tips of
  the needles of `camelshotnoise` does it, and the loop finds such vertices from the flattening
  instead of by hand. The probe, the threshold (20 in log units), the per-vertex rule, the drop to
  zero and the limits were chosen by measurement (crate README, "The feedback loop": on the
  dataset the clean meshes have `|δ|` at most 3.5 without cones and 7.6 after a solve, the needle
  tips of `camelshotnoise` `1.7·10⁴` without cones, and a cone at a tip hides it after a solve);
  the reasons are heuristic, not derived.
- **Why not formalized.** There is nothing to prove about a threshold. What is proved about the
  ingredients: the mismatch does not depend on the free constant of `u` or on the scale of the
  plane, and vanishes where the flattening realizes `u` up to those constants
  (`cone-feedback-mismatch`); the weights only decrease and every vertex is reweighted at most
  once (by reading: `is_reweighted`), so the loop ends after at most `max_rounds` solves.
- **Does correctness depend on it?** No, relative to the canonical problem **with the returned
  weights** (`FeedbackReport::weights`): every solve is a solve of a well-formed canonical problem
  (`wE ≥ 0`, at least one positive weight per closed component, guaranteed by the rule) and is
  certified as any other. The loop changes which problem is solved, not how. It is on by default in
  the pipeline (`ddg-uv`, Python, Blender); on every measured mesh other than the spiked ones it
  reweights nothing, so the result is that of the plain solve.

---

## Variational Surface Cutting

"Cuts paper" is Sharp, Crane, *Variational Surface Cutting*, ACM TOG 37(4), Article 156, 2018
(<https://doi.org/10.1145/3197517.3201356>). "CEA tutorial" is N. Sharp, *Céa's Method for
PDE-constrained shape optimization* (4/17/2018, updated 5/9/2019); "CEA derivations" is N. Sharp,
*Céa's Method Derivations for the Variational Shape Cuts Project* (May 8, 2019); both are the
author's notes. The contract of the implementation is `crates/ddg-variational-cuts/SPEC.md`
("SPEC"); its signs are those of `CLAUDE.md` (`Δ` positive semidefinite, `u` the true log scale
factor).

**What the machine-checked part covers for the cut flow, and what it does not.** There is no
canonical discrete problem here whose minimizer the code could certify: the flow is a
discretization of the gradient flow of a smooth energy. Machine-checked (`DdgProofs/Cuts/*.lean`,
ids `cut-*` in the README): every closed-form formula of one flow step is what SPEC says it is;
every linear system of a step has exactly one solution under the stated hypotheses; the discrete
adjoints are exact for the discrete objectives with respect to perturbations of the *source*;
the sign of the flow is a descent *given* the first-variation formula `dE = Σ g σ ℓ` and the
level-set relation `φ_t = −σ`; the local Hopf–Lax update is the minimum it is meant to be; and
the smooth formulas hold on three closed-form solutions (flat disk, spherical cap, shrinking
circle), where the shape derivative is computed directly by one-variable calculus; for the
default stages (SPEC §3.8–§3.10) and for `N` regions (§8): the well-posedness, conservation and
maximum principle of the systems on the cut, the convexity of the extension, the gradient norm,
and the combinatorial formulas of crossings, triple points and the projection. **Not**
machine-checked: the smooth results below (the four of the paper, the linearization behind the
roughening rate, the first variation of a network at a junction); that the discrete gradient density is a
consistent approximation of the shape derivative on a mesh (tested by finite differences and
refinement in `crates/ddg-variational-cuts/tests`, SPEC §6); that the discrete energy decreases
along the discrete flow; the accuracy claims of the flux recovery (SPEC D12 are measurements);
the marching order of the redistancing and its convergence; anything about rounding.

### [formula: cut-yamabe-equation] Yamabe equation with zero boundary values on the cut surface

- **Statement.** For a fixed cut `γ`, the Dirichlet energy `∫ |∇u|²` of the log conformal factor
  is minimized by a conformal flattening with `u = 0` on `∂M_γ`, and that `u` solves
  `Δu = −K` on `M_γ`, `u = 0` on `∂M_γ` (`Δ` positive semidefinite, `K` the Gaussian curvature).
- **Reference.** Cuts paper Sec. 3.1–3.2, Eq. 2–3; Springborn, Schröder, Pinkall 2008,
  *Conformal equivalence of triangle meshes*, App. E; Aubin 2013, *Some Nonlinear Problems in
  Riemannian Geometry*, Ch. 5.
- **Why not formalized.** Mathlib has no Gaussian curvature, Laplace–Beltrami operator or
  conformal-change formulas for Riemannian surfaces (see `smooth-yamabe-equation` above).
- **What depends on it.** It is the modelling step behind the linear system
  `u_B = 0`, `L_II u_I = −g_I` (SPEC §3.4). Proved instead: that system has exactly one solution
  (`[formula: cut-yamabe-system]`); `Σ_f Ω_f` is the total curvature (`cut-face-curvature-sum`);
  the first-order change of the discrete curvature under a conformal scaling is `L u`
  (`discrete-yamabe-linearization`, which fixes the sign); and the equation holds for the
  closed-form solutions on the disk and the spherical cap (`cut-disk-radial-solution`,
  `cut-sphere-cap-solution`: `diskU_ode`, `diskU_laplacian`, `capU_ode`).

### [formula: cut-shape-derivative] Shape derivatives of area and boundary integrals

- **Statement.** For a subdomain `Ω′` whose boundary moves with normal speed `σ`, and a function
  `ψ` that does not depend on `Ω′`: `D_σ ∫_{Ω′} dA = ∫_{∂Ω′} σ ds`,
  `D_σ ∫_{∂Ω′} ds = ∫_{∂Ω′} κ σ ds` (Eq. 6), `D_σ ∫_{Ω′} ψ dA = ∫_{∂Ω′} ψ σ ds` (Eq. 7),
  `D_σ ∫_{∂Ω′} ψ ds = ∫_{∂Ω′} (∂ψ/∂n + κψ) σ ds` (Eq. 8).
- **Reference.** Cuts paper Sec. 3.3.1, Eq. 6–8; CEA tutorial, section "Unconstrained Shape
  Derivative" (stated there "without derivation"); Allaire, Jouve, Toader 2004, *Structural
  optimization using sensitivity analysis and a level-set method*, J. Comput. Phys. 194.
- **Why not formalized.** Needs integration over moving domains of a surface (a transport
  theorem) and the first variation of arclength; Mathlib has neither for Riemannian surfaces.
- **What depends on it.** The term `β_r` of the region density (`Ddg.regionDensity`, Eq. 7 with
  `ψ = β`), and the length and penalty terms `∂w/∂n + wκ` of (★) in SPEC §1 (Eq. 8; SPEC D8 is
  the paper's omission of `∂η/∂n`), which the implicit step `(W_V + τ L_w) φ′ = W_V φ̃`
  discretizes. Machine-checked instances: `d/dρ (πρ²) = 2πρ` and `d/dρ ∫ ψ` for the radial
  integrands of the disk (`hasDerivAt_disk_regionEnergy`), the cap
  (`hasDerivAt_capEnergy_integral`), and the circle `d/dt (α·2πρ) = −(ακ)²·2πρ`
  (`hasDerivAt_circle_lengthEnergy`, `circle_levelSet_flow`).

### [formula: cut-cea-method] Céa's method and the shape derivatives of the distortion energies

- **Statement.** For `J(Ω′) = ∫_{Ω′} j(u)` with `Δu = f` in `Ω′`, `u = 0` on `∂Ω′` (Eq. 9–10),
  the shape derivative of `J` is the shape derivative of the Lagrangian
  `ℒ = ∫ j(u) + ∫ p(Δu − f) + ∫_{∂Ω′} λ u` at a critical point (Eq. 11), which gives
  `D_σ J = ∫_{∂Ω′} (j(u) + λ ∂u/∂n) σ ds` (Eq. 12). Table 1: Dirichlet `∫|∇u|²`: `Δp = −2Δu`,
  `λ = 0`, density `(∂u/∂n)²` (also Eq. 13); scaled Dirichlet `∫ a|∇u|²`: `Δp = −2Δ_a u`,
  `λ = −∂p/∂n − 2a ∂u/∂n`, density `−a(∂u/∂n)² − (∂u/∂n)(∂p/∂n)` (Eq. 18); pointwise `∫ b(u)`:
  `Δp = −b′(u)`, `λ = −∂p/∂n`, density `b(u) − (∂u/∂n)(∂p/∂n)`; `p = 0` on `∂Ω′` in all cases.
- **Reference.** Cuts paper Sec. 3.3, Eq. 9–13, Sec. 3.5, Eq. 17–19, Table 1 (the paper states the
  key step "which we will not justify here"); Céa 1986, *Conception optimale ou identification de
  formes, calcul rapide de la dérivée directionnelle de la fonction coût*, RAIRO Modél. Math.
  Anal. Numér. 20(3); CEA tutorial (the recipe); CEA derivations Sec. 1 (Dirichlet), Sec. 2
  (scaled Dirichlet), Sec. 3 (Hencky).
- **Why not formalized.** Shape calculus for PDE-constrained functionals (differentiability of
  `u` with respect to the domain, Sobolev spaces on surfaces, integration by parts on a moving
  domain); none of it is in Mathlib.
- **What depends on it.** The gradient densities `[formula: cut-gradient-dirichlet]`,
  `[formula: cut-gradient-scaled]`, `[formula: cut-gradient-pointwise]` and the adjoint systems
  `[formula: cut-adjoint-scaled]`, `[formula: cut-adjoint-pointwise]` are discretizations of
  Table 1; the direction of the flow depends on them. Proved instead: (1) the discrete adjoint
  is the exact first variation of the discrete objective with respect to the source
  (`scaledEnergy_perturb`, `hasDerivAt_pointwiseEnergy`); with `L_a = cL` it is `p = −2c u`,
  `q = −2c v` and the scaled density is `c` times the Dirichlet one (`scaledAdjoint_const`,
  `gradScaledD_const`); (2) on the flat disk and the spherical cap the derivative of the energy
  with respect to the radius, computed directly, **is** the Table 1 density integrated over the
  boundary, for the Dirichlet, the constant-weight scaled and the Hencky energy
  (`hasDerivAt_disk_dirichletEnergy`, `hasDerivAt_disk_scaledEnergy`,
  `hasDerivAt_disk_henckyEnergy` with the adjoint `diskP`, `hasDerivAt_capEnergy_integral`,
  `hasDerivAt_twoCaps` for two regions and their sign). Not proved: Table 1 for a general region,
  and that the discrete density converges to it (finite-difference tests, SPEC §6.3).

### [formula: cut-level-set-evolution] Evolution of the cut as a level set

- **Statement.** A function `φ` advected by a velocity field `X` satisfies `dφ/dt = −X·∇φ`
  (Eq. 20); for the signed distance function of `γ` with `∇φ = n` and the velocity `−σ* n` of the
  gradient flow this is `dφ/dt = σ*` on `γ` (Eq. 21); for curve shortening `X = −κ_φ n` it is
  `dφ/dt = −Δφ` with `Δφ = −κ_φ` (Eq. 22). Away from `γ` the speed is extended harmonically and
  `φ` is re-initialized to a signed distance function.
- **Reference.** Cuts paper Sec. 4.1–4.2, Eq. 20–22, Sec. 6, Fig. 9; Allaire, Jouve, Toader 2004;
  Osher, Sethian 1988, *Fronts propagating with curvature-dependent speed*, J. Comput. Phys. 79.
- **Why not formalized.** Level-set equations on surfaces (transport of a level set, mean
  curvature of level sets of a distance function) are not in Mathlib.
- **What depends on it.** The structure of one flow step (SPEC §4) and the relation `φ_t = −σ`
  between the change of `φ` and the normal speed of the cut, on which the sign of the flow rests.
  Proved instead: given that relation and the first variation `dE = Σ g σ ℓ`, the flow
  `φ_t = g` gives `dE = −Σ g² ℓ ≤ 0` and the opposite sign gives `+Σ g² ℓ`
  (`[formula: cut-speed]`, `firstVariation_levelSet`, `firstVariation_ascent`; SPEC D2); for a
  circle in the plane `φ = r − ρ(t)` with `ρ = √(ρ₀² − 2αt)` satisfies `φ_t = α(φ_rr + φ_r/r)` on
  the cut (`[formula: cut-curve-shortening-circle]`); for the two caps of a sphere the energy
  changes at the rate `(G_0 − G_1)·(length)` when the cut moves along the outward normal of
  region 0 (`hasDerivAt_twoCaps`). The harmonic extension and the implicit step are covered as
  linear algebra (`cut-harmonic-extension`, `cut-length-step`).

### [formula: cut-ripple-linearization] Growth rate of a ripple of the cut

- **Statement.** Linearizing the flow (★) of SPEC §1 about a straight cut between two half planes
  with constant normal derivatives, a ripple `δ(x) = ε cos(kx)` of the cut evolves by
  `δ_t = (ω k − w k²) δ`, `ω = 2 (Q_0 + Q_1)`, `Q_r` the part of the gradient density of region
  `r` that is quadratic in the normal derivatives. On a sphere of radius `R` with the Dirichlet
  energy the ripple `ε cos(mθ)` of the equator grows at the rate
  `(m − 1) (4/R³ − α_L (m + 1)/R²)`; the equator is stable iff `α_L > 4/(3R)`.
- **Reference.** SPEC §1 ("Roughening": the derivation, from the Dirichlet-to-Neumann map of a
  half plane resp. hemisphere, the second normal derivative of the base solution and the
  curvature of the ripple), §4, §6 item 5b (the rates were measured), D18. The Cuts paper states
  only the qualitative fact (Sec. 3.2.1: distortion alone is reduced by longer and longer cuts).
- **Why not formalized.** It is the linearization of a free-boundary problem: the derivative of
  the solution of a Poisson problem, and of its normal derivative, with respect to a perturbation
  of the domain. Nothing of this kind is in Mathlib.
- **What depends on it.** The choice of the length weight (`length_weight_for_ripples`,
  `[formula: cut-roughening-rate]`), the adaptive step (`time_scale/ω`), and the statements of
  SPEC about which ripples grow. Proved instead (`DdgProofs/Cuts/Ripple.lean`, `OnCut.lean`): the
  algebra of the rate (zero at the wavelength `ℓ` for `α_L = ωℓ/(2π)`, sign on either side,
  maximum, `ω ≥ 0` for the Dirichlet energy); the per-step factor of the discrete scheme
  (`rippleFactor_mem_Ioo`); the ingredients on closed forms: `e^{−ky} cos(kx)` is harmonic with
  normal derivative `k` times its boundary value (`halfPlane_mode_harmonic`,
  `halfPlane_mode_dtn`), the density difference `Q_0(1−kδ)² − Q_1(1+kδ)²` has the derivative
  `−ωk` (`hasDerivAt_ripple_density`), the mode `tan(ψ/2)^m` of the sphere has
  `sin ψ · h′ = m h` and the normal derivative `m` at the equator (`sphere_mode_ode`,
  `sphere_mode_dtn`), the base solution has `u″ = 1` there (`capU_second_derivative_equator`);
  and, from the cited formula, the stability threshold (`sphere_equator_stable_iff`,
  `sphere_unit_area_threshold`). Not proved: that these ingredients combine to the rate (the
  linearization), and the curvature `(m² − 1) ε/R²` of the rippled equator.

### [formula: cut-network-first-variation] First variation of the length of a network at a junction

- **Statement.** Minus the first variation of `∫ w ds` over a network of curves, for normal
  displacements of the curves and a displacement of a junction, is the curvature force of SPEC
  §3.10 on every curve and, at a junction, the sum of the unit tangents of the three curves
  leaving it (times `w`); for a constant weight it vanishes iff the curves meet at 120°.
- **Reference.** SPEC §8.5 ("Why"), D20; Mantegazza, Novaga, Tortorelli 2004, *Motion by
  curvature of planar networks*, Ann. Sc. Norm. Super. Pisa Cl. Sci. (5) 3. The Cuts paper
  (App. A) has no rule for junctions; the projection of the redistancing is from Losasso, Shinar,
  Selle, Fedkiw 2006, *Multiple Interacting Liquids*, ACM TOG 25(3), Sec. 3.1.
- **Why not formalized.** Calculus of variations on networks of curves on a surface.
- **What depends on it.** The rule `γ = 2π/3` at triple points of the default stages on networks
  (`[formula: cut-length-step-regions]`), which SPEC itself rates as no better than the paper's
  stages (D20: `FlowOptions::for_regions` is recommended). Proved instead
  (`DdgProofs/Cuts/Regions.lean`): the systems on the boundary of a region are well posed and
  keep the totals (`regionChain_*`); the two sides of a curve get opposite forces and rates
  (`curve_loads_opposite`, `onCut_rates_opposite`); the loads `2π/3 − β` of the three regions at
  a junction sum to zero and each vanishes at 120° (`junction_loads_sum`,
  `junction_load_balanced`). **Hypothesis of these, not derived here:** the interior angles of
  the regions at a cut point resp. a triple point sum to `2π` (the angles around a point inside
  a planar face).

### [formula: cut-fast-marching-convergence] The redistancing computes an approximate distance

- **Statement.** The values produced by the single-pass marching of SPEC §3.12 with the local
  Hopf–Lax update are a first-order approximation of the geodesic distance to the cut; on a mesh
  without obtuse angles the single pass is the unique solution of the fixed-point scheme
  `T_i = min over incident triangles`.
- **Reference.** Kimmel, Sethian 1998, *Computing geodesic paths on manifolds*, Proc. Natl. Acad.
  Sci. 95(15) (fast marching on triangulated surfaces; causality on acute triangles); Bornemann,
  Rasch 2006, *Finite-element discretization of static Hamilton–Jacobi equations based on a
  local variational principle*, Comput. Vis. Sci. 9(2) (the Hopf–Lax update, monotonicity and
  convergence); Cuts paper Sec. 6.4.
- **Why not formalized.** The marching order is an algorithm, not a formula, and convergence
  needs the theory of viscosity solutions of the eikonal equation.
- **What depends on it.** Nothing in the formulas: the redistancing only re-initializes `φ` away
  from the cut (cut-adjacent vertices keep their values). Proved: the local update is exactly
  `min_{λ∈[0,1]} λT_j + (1−λ)T_k + |x_i − (λx_j + (1−λ)x_k)|`
  (`[formula: cut-redistance]`, `hopfLaxTriangle_isLeast`, also from the side lengths alone),
  never exceeds the edge terms, and is at least `min(T_j, T_k)`. On obtuse triangles the single
  pass is the definition (SPEC D10).

## Navigating Intrinsic Triangulations

"Signpost paper" is Sharp, Soliman, Crane, *Navigating Intrinsic Triangulations*, ACM TOG 38(4),
Article 55, 2019 (<https://doi.org/10.1145/3306346.3322979>). The contract of the implementation
is `crates/ddg-intrinsic/SPEC.md` ("SPEC", decisions `N1`–`N23`). "B–S" is A. I. Bobenko,
B. A. Springborn, *A discrete Laplace–Beltrami operator for simplicial surfaces*, Discrete &
Computational Geometry 38(4), 740–756, 2007 (arXiv:math/0503219, which is what the Signpost paper
cites as "Bobenko and Springborn 2005").

**What the machine-checked part covers for intrinsic triangulations, and what it does not.**
Machine-checked (`DdgProofs/Intrinsic/*.lean`, ids `it-*` in the README): the closed forms of
single operations in exact real arithmetic. The layout of a triangle and of the two triangles on
an edge from their lengths is isometric and unique up to a rigid motion; the length of a flipped
edge, in three equivalent forms, the flip back, the triangle inequalities and the area of the new
triangles; strict convexity as crossing diagonals; the equivalence of the Delaunay condition and
a nonnegative cotan weight; that an edge which violates the condition can be flipped and
satisfies it afterwards; that a flip keeps the angle sums at its four vertices; the angle
normalization and the signpost update around a vertex; one step of a tracing query; barycentric
coordinates, the distances of a point to the corners and the lengths of an edge split; the
circumcentre in barycentric coordinates; the lengths and areas of a moved vertex's star; the
circumcentre and barycentre of a triangle times its area; the energy of Chen and Holst of one
star, its gradient (their Eq. 4.5, as a Fréchet derivative), that it is exactly quadratic, and
that the update of their Eq. 4.13 is a descent step (`DdgProofs/Intrinsic/Move.lean`). The theory of
one diamond (`DdgProofs/Delaunay/{Local,Rippa,Potential,Encroachment}.lean`): the equivalent local
Delaunay criteria (angles, cotangents, empty circumcircle of the layout), the flip of a
non-Delaunay edge with its corner-angle bookkeeping, Rippa's identity for one diamond
(`it-rippa-diamond`), two potentials that decrease strictly with each such flip
(`it-delaunay-potential`), and the encroachment criterion of Ruppert's rule
(`it-encroachment-angle`). The combinatorics of a triangulation and of the flip
(`DdgProofs/Delaunay/{Complex,Flip,Metric,FromTriMesh}.lean`, ids `it-complex-*`): a triangulated
surface as gluing data (the halfedge Δ-complex of `Complex`, with self-loops, multiple edges and
faces glued to themselves), its vertices, edges, boundary and Euler characteristic; that the edge
flip as `Complex::flip_edge` computes it keeps all of these, every vertex's class and cone angle,
a faithful vertex labelling and a metric (for a strictly convex diamond), and that two flips give
an isomorphic complex; Gauss–Bonnet on a closed one. **Not**
machine-checked: the eight results below; that the Rust arrays (`faces`, `twin`,
`halfedge_edge`, `edge_halfedge`, `vertex_halfedge`) are those of the gluing data (by reading, the
`it-complex-*` tags and `Complex::check`) and that the Rust code applies the formulas to the right
elements; a
tracing query as a whole (only one step is covered), edge traces and their checks (SPEC N7, N8),
the overlay (N12), point queries (N11), refinement and its termination (N1, N4, N14), the
simplicial export (N15); the shape limit of flips and every tolerance (N5, N6); anything about
rounding.

The sources below were checked on 2026-10-14 against primary or publisher-registry records; each
entry has a **Source check** line that says what was read, where, and what was not. "Read" means
the text of the work was read, not a summary of it. Bibliographic data (authors, title, venue,
volume, issue, pages, year, DOI) were taken from the Crossref record that the publisher
deposited for the DOI, and checked against the title page where the text was read. Section,
theorem and page numbers below refer to the copy that was read, named in the entry; where the
printed version may number differently, the entry says so. Where a statement could not be
checked in the work that is cited for it, the entry says that too and gives the work in which
the statement was found instead. The works that could be read in full: B–S (arXiv v3),
Indermitte et al. (open copy of the journal version, EPFL Infoscience), Chew 1993 (author's
copy of the proceedings paper), Shewchuk 1997 (CMU-CS-97-137) and 2002 (the Berkeley copy of
the 21 May 2001 manuscript of the journal article), Polthier–Schmies (the author's preprint),
and the Signpost paper. Not read: Rivin, Rippa, Fortune and the journal paper of Ruppert (paywalled,
no open copy found; of Ruppert only the first pages of the 1992 Berkeley technical report CSD-92-694,
an earlier version with another title, were read).

### [formula: it-delaunay-existence-uniqueness] The intrinsic Delaunay tessellation exists and is unique

- **Statement.** Every closed piecewise flat surface (a compact surface with a flat metric
  except at finitely many cone points) has a Delaunay tessellation with a given vertex set
  containing the cone points, and it is unique. Its cells are the images of the convex hulls of
  the vertices on the boundary of an immersed empty disk. A Delaunay triangulation, obtained from
  the tessellation by cutting every non-triangular cell by diagonals, is characterized by the
  local condition: a geodesic triangulation is Delaunay if and only if every interior edge is
  locally Delaunay, which is `θ_k + θ_l ≤ π` for the two angles opposite the edge
  (Signpost paper Eq. 4). Two Delaunay triangulations differ only in edges with `θ_k + θ_l = π`
  (vanishing cotan weight). Surfaces with boundary are reduced to the closed case by gluing
  flat caps to the boundary components.
- **Reference.** A. I. Bobenko, B. A. Springborn, *A discrete Laplace–Beltrami operator for
  simplicial surfaces*, Discrete & Computational Geometry 38(4), 740–756, 2007,
  <https://doi.org/10.1007/s00454-007-9006-1>, arXiv:math/0503219 (v3, 24 Feb 2006): Definitions
  2 and 3 (immersed empty disk, Delaunay tessellation), Proposition 4 (it is a tessellation;
  uniqueness is immediate from the definition, Sec. 2), Definition 8 and Lemma 9 (locally
  Delaunay if and only if the opposite angles sum to at most `π`), Proposition 10 (a
  triangulation is Delaunay if and only if all interior edges are locally Delaunay; the boundary
  by caps, in its proof). Existence by the flip algorithm: C. Indermitte, Th. M. Liebling,
  M. Troyanov, H. Clémençon, *Voronoi diagrams on piecewise flat surfaces and an application to
  biological growth*, Theoretical Computer Science 263(1–2), 263–274, 2001,
  <https://doi.org/10.1016/S0304-3975(00)00248-6>, Sec. 4, Theorem 1 (see the next entry).
  The statement was first made, without a proof of existence and uniqueness, by I. Rivin,
  *Euclidean structures on simplicial surfaces and hyperbolic volume*, Annals of Mathematics
  139(3), 553–580, 1994, <https://doi.org/10.2307/2118572>, Sec. 10 (as reported by B–S, Sec. 2).
  Signpost paper Sec. 4 (Eq. 4).
- **Source check.** Read, in B–S (arXiv v3): the statements and numbers above, and that the
  existence is "known" (B–S abstract) through the flip algorithm of Rivin and Indermitte et al.,
  while the uniqueness is theirs. The numbering is that of the arXiv version; the journal
  version was not compared. Read, in Indermitte et al.: Definitions 4 and 5 (legal edge, Delaunay
  triangulation by the local condition) and Theorem 1; they treat the boundaries of compact
  polyhedra (closed surfaces) and prove no uniqueness. Not read: Rivin (only the Crossref record:
  title, volume, issue, first page 553; the pages 553–580 and the section number are as B–S
  cite them). Corrected: the existence is no longer attributed to Indermitte et al. as a separate
  theorem, uniqueness to B–S only, Rivin as the first to state it; the boundary case
  (caps) is B–S's, not the Signpost paper's (the Signpost paper was not checked for it).
- **Proved, planar case (2026-11-05), existence only.** For a finite point set `S ⊂ ℝ²`:
  Delaunay's lemma, the planar form of B–S Proposition 10 (a triangulation of `S` that is locally
  Delaunay at every interior edge has empty circumcircles, `delaunay_lemma`; Delaunay ⟺ locally
  Delaunay, `isDelaunay_iff_isLocallyDelaunay`; empty circumdisk with the explicit circumcenter,
  `isDelaunay_iff_empty_circumdisk`), and existence of a Delaunay triangulation of `S` **given some
  triangulation of `S`** (`exists_delaunay_of_triangulation`; the existence of a first
  triangulation, e.g. by a sweep, is a hypothesis), `DdgProofs/Delaunay/PlanarLemma.lean`
  (`delaunay-lemma`). The planar local test equals layer 1's `θ_k + θ_l ≤ π` on the Euclidean
  lengths (`isLocallyDelaunay_dist2_iff`). Not proved: uniqueness up to ties, and the surface case.
- **Proved, planar case (2026-11-05), uniqueness up to ties.** Two Delaunay triangulations of a
  finite `S ⊂ ℝ²` differ only in triangles whose circumcircle passes through a fourth point of `S`
  (`exists_cocircular_of_not_mem`), so they are equal when no four points of `S` are cocircular
  (`NoFourCocircular`, required only for triples spanning a triangle; `delaunay_unique`); any two
  locally Delaunay triangulations are then equal whatever the flip order
  (`eq_of_isLocallyDelaunay`), and the Delaunay triangulation exists and is unique given a
  triangulation (`existsUnique_delaunay`), `DdgProofs/Delaunay/PlanarUniqueness.lean`
  (`planar-delaunay-uniqueness`; the lifting argument, de Berg et al. Thm. 9.6–9.7). Only the
  surface case of this entry stays cited.
- **Planar existence unconditional (2026-11-05).** Every finite `S ⊂ ℝ²` that is not collinear
  (three points with `sda ≠ 0`; necessary for `S ≠ ∅`) has a triangulation
  (`exists_triangulation`, `planar-triangulation-existence`; by induction, removing a
  lexicographically extreme point that keeps the rest non-collinear and coning it back to the
  visible boundary edges, `isTriangulation_insert`, `planar-triangulation-insertion`), hence a
  Delaunay triangulation (`exists_delaunay`, `planar-delaunay-existence`),
  `DdgProofs/Delaunay/Planar{Insertion,Existence}.lean`. The hypothesis of
  `exists_delaunay_of_triangulation` is discharged.
- **Proved, surface case: existence (2026-11-05, layer 4).** For every metric triangulation
  `s` (`MetricGluing`: closed or with boundary, any gluing), a Delaunay triangulation with the
  **original vertex set** is reachable by Delaunay flips (`MetricGluing.exists_delaunay`, and
  `flip_algorithm'` with the invariants χ, boundary, area, cone angles), `Delaunay/Surface/Finite.lean`
  (`it-delaunay-existence-uniqueness`; Indermitte et al. Thm. 1, B–S Thm. 15; termination as in
  `it-delaunay-flip-termination`). Delaunay ⟺ every interior edge locally Delaunay is the layer-3 definition
  (B–S Prop. 10 is the definition here; the empty-immersed-disk form is not modelled on a surface).
- **Still cited: the surface case of uniqueness.** Uniqueness of the Delaunay tessellation, and that
  two Delaunay triangulations of a surface differ only in edges with `θ_k + θ_l = π` (B–S, Prop. 4,
  Sec. 2); the empty immersed disks of B–S Definitions 2–3; the reduction of boundary surfaces to
  closed ones by caps. Layer 4's overlay route does not give them (they need empty immersed disks
  or connectivity of the flip graph of all geodesic triangulations with a given vertex set).
- **Why the rest is not formalized.** It is a global statement about all geodesic triangulations of a
  surface; neither piecewise flat surfaces nor immersed disks are in Mathlib.
- **What depends on it.** That the result of Delaunay flipping is *the* Delaunay triangulation up
  to ties, independent of the order of the queue (SPEC §5 item 4, N6: tested, the ties are
  counted). Proved without it: the local condition is equivalent to a nonnegative cotan weight
  (`[formula: it-cotan-weight]`, `cotanWeight_nonneg_iff`), and at a tie both diagonals satisfy
  it with equality (`flip_of_tie`). B–S Lemma 9 for one diamond (locally Delaunay ⟺ the
  apex `l` is not strictly inside the circumcircle of `i, j, k` in the unfolded diamond, and
  symmetrically) is proved in `DdgProofs/Delaunay/Local.lean`
  (`isLocallyDelaunay_iff_not_inside_circumcircle`, `inside_circumcircle_iff`, `on_circumcircle_iff`,
  from the power of the point, `2 A_jil (cot θ_k + cot θ_l)`).

### [formula: it-delaunay-flip-termination] Delaunay flipping terminates

- **Statement.** Starting from any geodesic triangulation of a compact piecewise flat surface,
  repeatedly flipping an interior edge that is not locally Delaunay ends after finitely many
  flips, in a Delaunay triangulation. Every such edge can be flipped (the two triangles form a
  convex quadrilateral, or the edge is already locally Delaunay; a topological obstruction is
  excluded). On a surface there is no bound on the number of flips in terms of the number of
  vertices alone. In the plane the number of flips is `O(n²)` (Signpost paper, citing Fortune).
- **Reference.** Termination: Indermitte et al. 2001 (as above), Sec. 4, Theorem 1: the sum of
  the areas of the circumcircles of the triangles decreases strictly with every flip and takes
  values in a discrete set. A second proof by Musin's harmonic index and Rippa's theorem, and the
  flipability of every non-Delaunay edge (Proposition 11): B–S (as above), Proposition 12 and
  Theorem 15, with the remark that Indermitte et al. missed the topological obstruction to
  flipping (B–S Sec. 2) and that Rivin's proof of termination is flawed. No bound on a surface:
  Indermitte et al., Sec. 4, after Theorem 1; Signpost paper Sec. 4.1. Planar bound: S. Fortune,
  *A note on Delaunay diagonal flips*, Pattern Recognition Letters 14(9), 723–726, 1993,
  <https://doi.org/10.1016/0167-8655(93)90142-Z>, as cited by the Signpost paper Sec. 4.1.
- **Source check.** Read: Indermitte et al. Theorem 1 and its proof (it refers to the Ph.D.
  thesis of H. Telley, 1989, Theorem 6.1, for the decrease of the circumcircle areas), Sec. 4;
  B–S Propositions 11 and 12, Theorem 15 (arXiv v3 numbering); Signpost Sec. 4.1. Corrected:
  the termination theorem is Indermitte et al.'s, with B–S giving the second proof and
  closing the gap (before: "B–S, whose abstract names the termination as proved there"); the
  pages of Fortune are added and the DOI is that of the Crossref record. Not checked: the content
  of Fortune's note. The `O(n²)` planar bound is quoted from the Signpost paper, which cites
  Fortune for it; Indermitte et al. (Sec. 2) cite for the same bound a technical report of
  H. Telley (EPFL, 1992), not Fortune. It is not used by the crate.
- **Formalized elsewhere (cross-reference only).** The planar case — Lawson flips of a plane
  triangulation of a fixed point set terminate in a Delaunay triangulation — is machine-checked in
  Coq: J.-F. Dufourd, Y. Bertot, *Formal study of plane Delaunay triangulation*, ITP 2010
  (arXiv:1007.3350; termination through a well-founded relation on a finite set, the decreasing
  quantity being the volume under the triangles lifted to the paraboloid). It is not a proof of the
  surface statement above (survey of 2026-10-31; the development's source could not be fetched).
- **Proved, planar case (2026-11-05).** For a finite point set `S ⊂ ℝ²` and triangulations of `S`
  in the point-set model (triangles with corners in `S` covering `conv S` with disjoint interiors,
  every point of `S` in a triangle a corner of it), `DdgProofs/Delaunay/Planar*.lean`: the Lawson
  flip of a strictly convex quadrilateral is again a triangulation of `S` (`lawson_flip`,
  `planar-lawson-flip`), an edge with the fourth corner strictly inside the circumcircle always
  has a strictly convex quadrilateral (`crossing_of_inCircumcircle`), the lifted volume
  `Σ area(t) Σ_{v∈t} |v|²` drops by exactly `sda·inCircle/(2|sda|) > 0` at every such flip
  (`liftPotential_flip`, `planar-flip-potential`), there are at most `2^(|S| choose 3)`
  triangulations, so the flip relation is well founded and every flip sequence has at most that
  many steps (`wellFounded_flipStep`, `flip_sequence_length_le`, `planar-flip-termination`), and it
  ends in a locally Delaunay — by Delaunay's lemma (`delaunay_lemma`, `delaunay-lemma`) a
  Delaunay — triangulation (`exists_flips_to_delaunay`). Planar only; the surface statement and
  its finiteness hypothesis below were unchanged by it (the surface finiteness is now proved, layer 4).
- **Proved, planar `O(n²)` bound (2026-11-05).** Every sequence of Delaunay flips from a
  triangulation of `S` has at most `C(|S|, 2)` flips (`flip_count_le`), because an edge flipped
  away never reappears (`flipped_edge_not_mem`): the lifted surface goes down with every flip
  (`liftBelow_flip`) and lies strictly below the lifted segment of the flipped edge at the crossing
  of the diagonals, while on an edge it equals the lifted segment (`liftPlane_eq_of_mem_edge`),
  `DdgProofs/Delaunay/PlanarFlipCount.lean` (`planar-flip-count`; Edelsbrunner, *Geometry and
  Topology for Mesh Generation*, Sec. 1.3). Fortune's note was not used.
- **Proved (trimmed, 2026-11-05, layer 3).** Everything but the finiteness, on the gluing data of
  `DdgProofs/Delaunay/Complex.lean` (self-loops, multiple edges, faces glued to themselves):
  `DdgProofs/Delaunay/Termination.lean` (`it-delaunay-self-glued`, `it-delaunay-potential`,
  `it-delaunay-flip-algorithm`) and `DdgProofs/Delaunay/FlipQueue.lean` (`it-delaunay-flip-loop`).
  Every non-Delaunay interior edge can be flipped (an edge with one face on both sides is always
  Delaunay, strictly: `IsMetric.locallyDelaunayAt_of_same_face`; B–S's topological obstruction);
  the global lifting energy and the harmonic index decrease strictly with every Delaunay flip
  (`liftPotential_flip_lt`, `harmonicPotential_flip_lt`), so no state repeats; if the states
  reachable by Delaunay flips are finitely many, every flip sequence is finite
  (`MetricGluing.acc_of_finite`, at most `|reach| − 1` flips), and its end is Delaunay with the
  χ, boundary, area and cone angles of the start (`MetricGluing.flip_algorithm`). The literature's
  hypothesis — for every `L`, finitely many triangulations of the flip class with edges `≤ L`, up
  to isomorphism (`MetricGluing.BoundedLengthFinite`) — implies that finiteness
  (`reach_finite_of_boundedLengthFinite`: edges stay `≤ √(H·A)` because `H` decreases and the
  area is constant, B–S's argument). The Rust queue loop with its tolerance is proved to make only
  Delaunay flips, to terminate under the same hypothesis in any order, and, when it exits without
  the cap, to leave every interior edge passing its test (`MetricGluing.loop_result`).
- **Proved, layer 4 (2026-11-05): the finiteness, and so unconditional termination on a surface.**
  `MetricGluing.boundedLengthFinite` (the former hypothesis `BoundedLengthFinite`) is a theorem for
  every `MetricGluing` (any closed surface or surface with boundary, with self-loops, multiple edges
  and self-glued faces), `DdgProofs/Delaunay/Surface/*.lean` (`Finite.lean`; plan and route in
  `Delaunay/LAYER4-PLAN.md`). The cone surface and its geodesics are not modelled: every
  triangulation reached by flips is an *overlay* on the original one (faces laid out in the plane
  together with the faces of `s` that cover them; a convex flip keeps this, `Overlay.flip`), every
  edge of it is a straight chain through at most `crossBound s L` faces of `s` (`Overlay.exists_word_le`),
  and there are finitely many such crossing words (`Word.finite_length_le`). Consequently
  `reach_finite`, `acc` (**every sequence of Delaunay flips is finite**, from every state),
  `wellFounded_delaunayFlip`, `not_exists_infinite_flips'`, `flips_length_le'` (at most
  `|reach s| − 1` flips), `exists_delaunay`, `flip_algorithm'` and, for the Rust loop in any
  order, `loop_acc'` hold without hypothesis. Only the remark that no bound in `|V|` alone exists
  (Statement above) stays unproved; it is not a theorem used anywhere. Nothing else of this
  entry is cited, apart from the sources and the planar `O(n²)` bound of Fortune's note.
- **What depends on it.** That `Delaunay flipping` returns at all with `FlipLimit::Unlimited`;
  the crate has a limit on the number of flips and reports when it is reached (SPEC §4, N6; never
  reached in the recorded runs). Earlier one-diamond results (`diamondConvex_of_not_delaunay`,
  `delaunay_of_same_triangle`, `flip_of_not_delaunay`, `liftEnergy_flip`, `harmonicIndex_flip_lt`)
  are the steps of the proof above.

### [formula: it-rippa-theorem] Rippa's theorem: Delaunay flips do not increase the Dirichlet energy

- **Statement.** For data `f` at the vertices, let `f_T` be the piecewise linear interpolant on a
  triangulation `T`. If `T₂` arises from `T₁` by flipping an edge `e` and `e` is locally
  Delaunay in `T₂`, then the Dirichlet energy satisfies `E(f_{T₁}) ≥ E(f_{T₂})`, with equality only
  if `f_{T₁} = f_{T₂}` or `e` was already locally Delaunay in `T₁`. Hence the minimum of the
  energy over all triangulations is attained on the Delaunay triangulations. This holds in the
  plane (Rippa) and, with the same proof, for geodesic triangulations of a piecewise flat surface
  (B–S), so the cotan Laplacian of a Delaunay triangulation depends on the surface and its
  vertices only.
- **Reference.** S. Rippa, *Minimal roughness property of the Delaunay triangulation*, Computer
  Aided Geometric Design 7(6), 489–497, 1990, <https://doi.org/10.1016/0167-8396(90)90011-F>;
  for piecewise flat surfaces B–S (as above), the statement "Rippa's Theorem" in Sec. 2, between
  Definition 13 and Proposition 14, and Sec. 3 (Definition 16, Proposition 17) for the Laplacian. The Signpost paper mentions
  the theorem in Sec. 5.5 (Globally Optimal Direction Fields), as a contrast to vector fields.
- **Source check.** Read: the statement in B–S (arXiv v3), including their remark that Rippa's
  proof holds without change for piecewise flat surfaces, and the Signpost paper's use. Not read:
  Rippa's paper (paywalled): its own wording, theorem number and the planar statement are taken
  from B–S's quotation; bibliographic data from Crossref. Corrected: the statement follows
  B–S's form (the condition is on the edge after the flip, with the case of equality) instead of
  "an edge that violates the condition".
- **Proved: the one-diamond part** (`DdgProofs/Delaunay/Rippa.lean`, `it-rippa-diamond`). For a
  strictly convex diamond the difference of the cotan energies of the two triangulations is
  exactly `−(cot θ_k + cot θ_l)/(2 A_ilk A_jkl)·ℓ(u)²`, with `ℓ` the affine dependence of the four
  vertices (`rippa_identity`); hence the flip of an edge with `θ_k + θ_l ≥ π` does not increase
  the energy for any `u` (`rippa_le`), strictly if `θ_k + θ_l > π` and `ℓ(u) ≠ 0` (`rippa_lt`),
  with equality at a tie (`rippa_eq_of_tie`), and for a convex diamond this holds for all `u`
  only if `θ_k + θ_l ≥ π` (`forall_rippa_le_iff`). `cotEnergy` is the Dirichlet energy of the
  linear interpolant on each triangle by `cotan_formula`.
- **Why not formalized (the rest).** The global statement (minimum over all triangulations) needs
  the previous entry and a model of all triangulations of a surface with a given vertex set.
- **What depends on it.** Only the remark of SPEC N6 that two results of flipping that differ at
  ties have the same Dirichlet energy for every vertex function (tested). Proved here: the
  off-diagonal entry of the cotan Laplacian at an edge with two faces is minus half the sum of
  the two opposite cotangents (`offDiag_eq_neg_half_cotanWeight`) and is nonpositive exactly
  when that sum is nonnegative (`offDiag_nonpos_iff`), which for two opposite angles in
  `(0, π)` is the Delaunay condition (`cotanWeight_nonneg_iff`).

### [formula: it-chew-second-algorithm] Chew's second algorithm in the plane

- **Statement.** In the plane, Delaunay refinement that inserts the circumcentre of a triangle
  with an angle below the bound into the constrained Delaunay triangulation terminates. Chew's
  algorithm for a bound of `30°` assumes input angles of at least `60°` (angles down to `30°` with
  a modification, smaller ones are cut off); a circumcentre that would lie beyond a boundary
  segment is not inserted, the segment is split instead (in halves, in thirds if its length is
  in a short range) and the free vertices closer than the new edge length to the new vertex are
  deleted. The argument: an angle below `30°` is opposite a side shorter than the circumradius,
  so the circumradius of a failing triangle is larger than an existing edge and no edge shorter
  than the shortest one present is created. Shewchuk's reanalysis of the simpler form of the
  rule (split the encroached segment at its midpoint and delete the free vertices inside its
  diametral circle) terminates for angle bounds up to `arcsin(1/√5) ≈ 26.56°` (circumradius to
  shortest edge ratio `B ≥ √5/2`) when input angles are at least `60°`, and then also gives good
  grading and size optimality.
- **Reference.** L. P. Chew, *Guaranteed-quality mesh generation for curved surfaces*, Proc. 9th
  Annual Symposium on Computational Geometry (SCG '93), ACM, 1993, 274–280,
  <https://doi.org/10.1145/160985.161150>: Sec. 2 (Algorithm A, with the assumption on input
  angles on p. 275, step 6, Lemma 2 on p. 276, Theorem 3 on p. 277: "Algorithm A halts" and all
  angles at least `30°`). J. R. Shewchuk, *Delaunay refinement algorithms for triangular mesh
  generation*, Computational Geometry 22(1–3), 21–74, 2002,
  <https://doi.org/10.1016/S0925-7721(01)00047-5>: Sec. 4 (Chew's second algorithm), 4.1 (the
  description above), 4.2 (termination for `B ≥ √5/2`, Lemma 5, `26.56°`), 4.3 (grading), 7.3
  (Chew's own bound `B = 1`, `30°`, by range-restricted segment splitting). J. R. Shewchuk,
  *Delaunay Refinement Mesh Generation*, Ph.D. thesis, Carnegie Mellon University, CMU-CS-97-137,
  1997: Sec. 3.4 and 3.5.2 (the same), Sec. 5.3.2, pp. 138–139: "All of the proofs in Chapter 3
  can be made to work with geodesic distances" for planar domains with holes and for flat
  surfaces meeting along segments in three dimensions; this is the passage that the Signpost
  paper (Sec. 4.2) cites as only suggesting that the lemmas carry over.
- **Source check.** Read, as listed in the preamble: Chew 1993 Sec. 2 and the start of Sec. 3
  (which rejects geodesic circles on curved surfaces and defines the circumcircle through
  spheres, so the paper says nothing about cone surfaces); Shewchuk 2002 (the 2001 manuscript:
  numbers are those of its sections) and 1997; Signpost Sec. 4.2. Corrected: `30°` is Chew's
  bound and `26.5°` is Shewchuk's (before: Chew's algorithm "for angle bounds up to about
  `26.5°`" with `30°` only without boundary); the boundary rule is Chew's step 6, and the
  diametral-circle form is Shewchuk's description (Sec. 4.1) of it; the assumption on input
  angles is added. The journal version of Shewchuk 2002 was not compared with the manuscript. In
  the Signpost paper the claim for planar triangulations is cited to Chew 1987 (*Constrained
  Delaunay triangulations*, SCG '87), which was not read.
- **Why not formalized.** A termination proof for an algorithm on planar triangulations by a
  packing argument; no planar Delaunay triangulations in Mathlib.
- **What depends on it.** Nothing that the crate claims: for intrinsic triangulations SPEC N1
  says what is *argued* (closed surfaces, every angle sum at least `60°`, bounds up to `30°`) and
  what is only *observed* (with boundary; with vertices of small angle sum, N4), and the limit on
  insertions is reported when reached. Proved here, of the planar lemma's ingredients: the
  circumcentre in barycentric coordinates is equidistant from the three corners at the distance
  `a₀a₁a₂/4A` (`circumcenter_equidistant`), it lies beyond a side exactly when the opposite
  angle is obtuse (`circumcenterBary_neg_iff`), and an angle below `30°` is opposite a side
  shorter than the circumradius (`circumradius_gt_of_angle_lt`; also, from lengths only,
  `DdgProofs/Delaunay/Encroachment.lean`: `circumradius_eq_div_sin`,
  `lt_circumradius_of_lt_pi_div_six`).

### [formula: it-ruppert-encroachment] Ruppert's boundary rule in the plane

- **Statement.** In the plane, Delaunay refinement that splits a boundary segment whenever a
  vertex lies on or inside its diametral circle and is visible from it (the segment is
  *encroached*), and inserts a circumcentre only if it encroaches on no segment (otherwise the
  circumcentre is rejected and the encroached segments are split), terminates for angle bounds up
  to `arcsin(1/(2√2)) ≈ 20.7°` (circumradius to shortest edge ratio `B ≥ √2`) when no two
  input segments meet at an angle below `60°` (Shewchuk's analysis); it then creates no edge
  shorter than the minimum local feature size of the input. For smaller input angles, splitting
  the segments at an input vertex on concentric circular shells (radii that are powers of
  two) prevents the two segments from encroaching on each other for ever. (The criterion that a
  vertex inside the diametral circle of a segment subtends an obtuse angle is proved, see below.)
- **Reference.** J. Ruppert, *A Delaunay refinement algorithm for quality 2-dimensional mesh
  generation*, Journal of Algorithms 18(3), 548–585, 1995,
  <https://doi.org/10.1006/jagm.1995.1021>. For the statement as given: J. R. Shewchuk, 2002
  (as above), Sec. 3.1 (the rules), Theorem 4 in Sec. 3.3 with the remark after its proof
  (`B = √2` gives `20.7°`), Sec. 6.2 (the concentric shells, attributed to Ruppert); J. R. Shewchuk,
  1997 (as above), Sec. 3.3, and Sec. 5.3.2 (encroached segments are diagnosed by an obtuse
  angle opposite them).
- **Source check.** Read: Shewchuk 2002 and 1997 at the places given. Not read: Ruppert's
  journal paper (paywalled). The first pages of its earlier version, the 1992 technical report
  (UC Berkeley CSD-92-694, titled *A New and Simple Algorithm for Quality 2-Dimensional Mesh
  Generation*), state only the aspect ratio and size guarantees of the algorithm, so the
  `20.7°` and the condition on input angles were not seen in Ruppert's own text. Corrected:
  before, the entry required no input angle below `90°` ("below `60°` in Shewchuk's
  analysis"); only the `60°` of Shewchuk's Theorem 4 was found in a primary source read here (a
  secondary source, Pav and Walkington, reports that Ruppert's analysis needs non-acute input
  angles, which was not checked against Ruppert), so the `90°` is dropped. The encroachment
  condition now says "on or inside".
- **Why not formalized.** As the previous entry.
- **What depends on it.** SPEC N1 uses the rule on surfaces with boundary and states that its end
  at `25°` and `30°` is observed, not proved, as in the plane. Proved here
  (`DdgProofs/Delaunay/Encroachment.lean`, `it-encroachment-angle`, the test of `is_encroached`):
  a point is strictly inside the diametral circle of a segment iff the angle it subtends is
  obtuse, on it iff the angle is right (`dist_midpoint_lt_iff`, `dist_midpoint_eq_iff`), and in
  lengths, the corner opposite `s` is obtuse iff `x² + y² < s²` (`pi_div_two_lt_cornerAngle_iff`,
  `encroaches_iff_cornerAngle`); the criterion is also the sign of a coordinate of the
  circumcentre (`circumcenterBary_neg_iff`). Not proved: termination of the rule.

### [formula: it-straightest-geodesic] Tracing by unfolding follows the straightest geodesic

- **Statement.** On a polyhedral surface, a *straightest geodesic* is a curve whose left and
  right curve angles `θ_l`, `θ_r` (with `θ_l + θ_r` the total angle of the point) are equal at
  every point. Given a point and a polyhedral tangent vector there is a unique straightest
  geodesic with that start point and direction, and it extends to the boundary of the surface.
  A geodesic that contains no vertex is both straightest and locally shortest, and it is the
  straight line in the unfolding of the faces along its way; through a vertex the two notions
  differ (a straightest geodesic through a spherical vertex is not locally shortest).
- **Reference.** K. Polthier, M. Schmies, *Straightest geodesics on polyhedral surfaces*, in
  H.-C. Hege, K. Polthier (eds.), Mathematical Visualization, Springer, 1998, 135–150,
  <https://doi.org/10.1007/978-3-662-03567-2_11>: Definition 8 (straightest geodesic), Theorem 10
  (discrete initial value problem), Lemma 11 (1) and its proof (unfolding the faces along the
  curve), 11 (2)–(3). Signpost paper Sec. 3.2.2, which names the tracing query as the discrete
  exponential map producing a straightest polyhedral geodesic in this sense.
- **Source check.** Read: the author's preprint (18 pages, Polthier's page at FU Berlin,
  `straightest_preprint.pdf`); the numbers above are those of the preprint, the printed chapter
  (pages and editors from Crossref: H.-C. Hege, K. Polthier; Springer Berlin Heidelberg, 1998) was
  not compared. The Signpost sentence was read. Corrected: the statement now says in terms of
  the curve angles what "straight" means and gives the unfolding as Lemma 11's proof; before,
  "unique through vertices" and "equal angle sums" were stated without a number.
- **Why not formalized.** It is a statement about curves on a polyhedral surface, which is not
  formalized here; a path through many triangles is a composition of the steps below with the
  bookkeeping of the mesh.
- **What depends on it.** That a tracing query evaluates the exponential map of the surface.
  Paths through input vertices are not continued by the crate (SPEC N7: the trace stops there);
  through inserted vertices, which are flat, it goes straight (N11). Proved here, for one step:
  the unfolded triangle has the given lengths and lies on the far side of the entry edge
  (`unfoldApex_isometric`), the path leaves through exactly one of the two other edges unless it
  hits the opposite vertex (`exit_cases`), the crossing is on that edge
  (`axisCrossingT_mem_Icc`) and lies ahead (`exit_ahead_above`, `exit_ahead_below`).

### [formula: it-odt-delaunay-optimality] The Delaunay connectivity minimizes the interpolation error; the boundary rule

- **Statement.** For a point set `V` of the plane (convex hull `Ω`) and `u(x) = |x|²`, the error
  of the piecewise linear interpolant `u_I` on a triangulation `T`, `E(T) = ∫_Ω |u_I − u| dx`, is
  smallest for the Delaunay triangulation of `V` among all triangulations of `V` (Chen–Holst
  Theorem 2.3, proved in Chen–Xu; first in the plane by Rippa for `L^q`, `1 ≤ q < ∞`, and by
  D'Azevedo–Simpson for `q = ∞`). An *optimal Delaunay triangulation* minimizes `E` over the
  positions of the interior vertices as well. For *boundary elements* (simplices with at least
  one vertex on the boundary) Chen–Holst replace the circumcentre by the barycentre (Sec. 4.2.1,
  "For boundary elements, use the barycenter"), a remedy against stretched triangles at the
  boundary, not a gradient. The Signpost paper (Sec. 4.3) applies Eq. 4.13 to inserted vertices
  of an intrinsic triangulation, averaging the *vectors* to the circumcentres in the tangent
  plane of the vertex, followed by Delaunay flipping.
- **Reference.** L. Chen, M. Holst, *Efficient mesh optimization schemes based on Optimal
  Delaunay Triangulations*, Computer Methods in Applied Mechanics and Engineering 200(9–12),
  967–984, 2011, <https://doi.org/10.1016/j.cma.2010.11.007>: Theorem 2.3, Lemmas 4.1–4.3
  (Eq. 4.2–4.8), Eq. 4.13 and the paragraphs after it, Sec. 4.2.1 (numbering of the preprint,
  see below). L. Chen, J. Xu, *Optimal Delaunay triangulations*, Journal of Computational
  Mathematics 22(2), 299–308, 2004, <https://doi.org/10.4208/jcm.v22.n2.p299>. P. Alliez,
  D. Cohen-Steiner, M. Yvinec, M. Desbrun, *Variational tetrahedral meshing*, ACM Transactions on
  Graphics 24(3), 617–625, 2005 (SIGGRAPH 2005), <https://doi.org/10.1145/1073204.1073238>.
  Signpost paper Sec. 4.3.
- **Source check.** Read: the preprint of Chen–Holst (`../NavigatingIntrinsicTriangulations/cited/
  chen-holst-2011/ChHo10a.pdf`, PDF dated 2 Nov 2010, title "Efficient Mesh Optimization Schemes
  based on Optimal Delaunay Triangulations"): Theorem 2.3 and Remark 2.4, Sec. 4.2.1 with
  Lemmas 4.1–4.3, the proof of Eq. 4.5, Eq. 4.13, the line search `α = 1, 1/2, 1/4, …` and the
  boundary rule; the attribution of Eq. 4.3 to Chen's 2004 Meshing Roundtable paper and of
  Eq. 4.5 to Alliez et al. is the preprint's. The printed version was not compared, so section
  and equation numbers may differ there. Bibliographic data of all three works from their Crossref
  records (authors, titles, venues, volume 200, issue 9–12, pages 967–984, February 2011;
  volume 22, issue 2, pages 299–308, 2004; volume 24, issue 3, pages 617–625, July 2005), and the
  volume and pages of Chen–Xu and Alliez et al. also from the preprint's reference list ([2],
  [24]). Not read: Chen–Xu and Alliez et al. (only their Crossref records); the statements
  attributed to them are as Chen–Holst report them. Chen–Holst do **not** claim that Eq. 4.13
  decreases the energy for every `α`; the line search guards validity of the mesh, not descent.
- **Proved instead** (formerly part of this entry, id `it-odt-energy`;
  `DdgProofs/TriangleIntegral/Tiling.lean`, with `it-odt-interpolation-error` and `it-odt-update`):
  Lemma 4.1 for a triangulation, as a Lebesgue integral over any tiling of `Ω` by the triangles
  (`IsTiling`: union `Ω`, disjoint interiors; boundaries of triangles are null sets),
  `E = Σ_T (|T|/3) Σ_i |x_i|² − ∫_Ω |x|² = (1/3) Σ_i |x_i|² |ω_i| − ∫_Ω |x|²
  = Σ_T (|T|/12)(ℓ₀₁² + ℓ₁₂² + ℓ₂₀²)`; for a fan around `x` plus fixed triangles,
  `E(x) = S(x) + C − ∫_Ω |y|²` with the star energy `S = odtEnergy` of `Intrinsic/Move.lean`, whose
  derivative (Eq. 4.5), Hessian (Lemma 4.3) and the exact change under the step of Eq. 4.13 are
  proved there; so the step with `0 ≤ α ≤ 2` lowers `E` itself and `α = 1` minimizes it, whenever
  the moved fan still tiles the same region (`odt_error_step_le`). The tiling is a hypothesis: that
  the moved fan tiles the same star region is not proved here (the line search of N21 keeps every
  triangle of the star positively oriented, which the theorems also assume).
- **Why not formalized.** Theorem 2.3 is a statement about all triangulations of a point set (the
  lifting of `V` to the paraboloid and the lower convex hull, as for
  `it-delaunay-existence-uniqueness`); triangulations of a point set are not modelled in this
  tree and Mathlib has no Delaunay triangulation. The boundary rule is a heuristic, not a theorem.
- **What depends on it.** The update of the optimal Delaunay triangulation (SPEC §4, N21) and the
  barycentres for triangles at the boundary. That flipping to Delaunay after a move does not
  increase `E` is Theorem 2.3 in the plane; for an intrinsic triangulation it is the intrinsic
  analogue (`it-rippa-theorem`, Rippa's theorem being the `L^q` case). Not claimed: that the
  intrinsic iteration converges, or that a step shortened to stay in the star (N19, N20) keeps
  the descent beyond what `odt_error_step_le` states for a planar star.

### [formula: it-odt-target-length-splits] Splitting edges longer than a target length between optimization steps

- **Statement.** Vertex optimization alone cannot add vertices; the input vertices of an
  intrinsic triangulation cannot move. So between iterations of the optimal Delaunay update the
  Signpost paper inserts new vertices by splitting edges longer than a user-defined target
  length, "à la Tournois et al. [2009]" (Sec. 4.3). Tournois et al. interleave Delaunay
  refinement (Steiner points inserted while an element violates the criteria; an edge is *badly
  sized* if it is longer than the sizing function at its midpoint) with passes of ODT
  optimization, optimizing only after a refinement batch inserted points.
- **Reference.** J. Tournois, C. Wormser, P. Alliez, M. Desbrun, *Interleaving Delaunay
  refinement and optimization for practical isotropic tetrahedron mesh generation*, ACM
  Transactions on Graphics 28(3), Article 75, 2009 (SIGGRAPH 2009),
  <https://doi.org/10.1145/1531326.1531381>; research report INRIA RR-6826, 2009,
  <https://inria.hal.science/inria-00359288>. Signpost paper Sec. 4.3.
- **Source check.** Read: the INRIA report RR-6826 (HAL inria-00359288v1, submitted 6 Feb 2009):
  Sec. 2 (Algorithm 1: refine, then optimize, while refinement inserted points), Sec. 2.2 and 2.4
  (the sizing criterion: an edge longer than the sizing function at its midpoint is bad),
  Sec. 2.5 (ODT optimization, citing Chen 2004 and Alliez et al. 2005). Bibliographic data from
  the Crossref record (authors, title, volume 28, issue 3, July 2009; Crossref gives pages 1–9
  and no article number); the article number 75 from the Signpost paper's reference list
  ("28, 3 (2009), 75"). Not read: the journal version. The report treats tetrahedral meshing in
  a domain of `ℝ³` with Delaunay refinement of bad elements, not the splitting of edges at their
  midpoints; the Signpost paper's rule is its own adaptation, for which Tournois et al. are the
  model, not the source of a theorem.
- **Why not formalized.** It is a heuristic of an algorithm, with no statement to prove; what is
  checked for a split (lengths, signposts) is `it-point-distances`.
- **What depends on it.** The schedule of the optimal Delaunay triangulation (SPEC §4, N21, N22):
  which edges are split, where, and how often. Nothing of correctness depends on it; the
  vertices it adds are inserted vertices like any other (N10).

## Navigating Intrinsic Triangulations: applications 5.1 to 5.3

The applications of the Signpost paper's Sec. 5.1 (Steiner trees), 5.2 (finite elements) and 5.3
(geodesic distance by the heat method) are in `crates/ddg-intrinsic-apps` (contract: its `SPEC.md`,
decisions `N30`–`N39`). Machine-checked (`DdgProofs/IntrinsicApps/HeatFem.lean`): the closed forms
of one triangle, `app-heat-gradient`, `app-heat-divergence`, `app-fem-galerkin-mass` (row sums and
the bounds against the lumped mass; that its entries are the integrals `∫ φ_i φ_j` is proved in
`DdgProofs/TriangleIntegral/Fem.lean`) and `app-fem-midpoint-subdivision`. The stiffness matrix is
`cotan-fem-stiffness`, proved as the matrix of the integrals `∫ ∇φ_i·∇φ_j` in the same file; its
convergence is `cotan-fem-convergence` above. Cited below: the Fermat point of a triangle. Not covered: the Steiner
heuristic as a whole (Smith, Lee, Liebman 1981, a heuristic with no optimality statement), the
convergence of the heat method (Crane, Weischedel, Wardetzky 2013, Sec. 3.2 and their
Appendix A), the exact window propagation of `exact.rs` (Mitchell, Mount, Papadimitriou 1987;
Surazhsky et al. 2005; the filter of Xin and Wang 2009, argued in its module documentation), the
Lanczos iteration, and every tolerance.

### [formula: app-steiner-fermat-point] The Fermat point of a triangle and the Steiner tree of its corners

- **Statement.** In a triangle whose angles are all below `120°` (`2π/3`) the point `F` that
  minimizes the sum of the distances to the three corners (the Fermat or Torricelli point, the
  first isogonic centre) lies inside the triangle and sees every side under `120°`; the shortest
  network joining the three corners is the three segments from `F` (the Steiner minimal tree of
  three points; if an angle is `120°` or more, it is the two sides at that corner). Its
  trilinear coordinates are `csc(A + π/3) : csc(B + π/3) : csc(C + π/3)`, so its barycentric
  coordinates are `a csc(A + π/3) : b csc(B + π/3) : c csc(C + π/3)` (`a` the side opposite the
  angle `A`), and the length of the tree is `S` with `S² = (a² + b² + c²)/2 + 2√3·Area`.
- **Reference.** E. N. Gilbert, H. O. Pollak, *Steiner minimal trees*, SIAM Journal on Applied
  Mathematics 16(1), 1–29, 1968, <https://doi.org/10.1137/0116001> (the three-point problem and
  the `120°` condition). C. Kimberling, *Encyclopedia of Triangle Centers*, X(13) (Fermat point;
  trilinears `csc(A + π/3)`), <https://faculty.evansville.edu/ck6/encyclopedia/ETC.html>.
  Signpost paper Sec. 5.1, which prints "barycentric coordinates `b_i := csc(θ_i + π/3)`".
- **Source check.** Not read: Gilbert–Pollak (paywalled) and the Encyclopedia entry; the
  statements are the classical ones and are checked numerically instead
  (`crates/ddg-intrinsic-apps/tests/steiner.rs`, `the_fermat_point_joins_the_corners_at_120_degrees`:
  the three directions from the computed point meet at `120°` to `1e-12` and the length agrees
  with the formula for `S²` to `1e-12`, on four triangles). The paper's `csc(θ_i + π/3)` used as
  barycentric coordinates gives another point (they are the trilinear coordinates); SPEC N30.
- **Why not formalized.** The minimality of the network is a statement about all networks
  joining three points; the coordinates need the trigonometry of the triangle's angles, for
  which nothing in `lean/proofs` is set up.
- **What depends on it.** The length and the place of every Steiner point of `steiner_tree`
  (SPEC N30). Nothing else: the tree is checked as a tree and the cut as a disk whatever the
  points are.

---

## Navigating Intrinsic Triangulations, Sec. 5.4 and 5.5: adaptive refinement and tangent fields

The contract is `crates/ddg-intrinsic-apps/SPEC.md` (decisions N40–N49). "VHM" is N. Sharp,
Y. Soliman, K. Crane, *The Vector Heat Method*, ACM TOG 38(3), Article 24, 2019
(<https://doi.org/10.1145/3243651>); "GODF" is F. Knöppel, K. Crane, U. Pinkall, P. Schröder,
*Globally Optimal Direction Fields*, ACM TOG 32(4), Article 59, 2013
(<https://doi.org/10.1145/2461912.2462005>); "MNS" is P. Morin, R. H. Nochetto, K. G. Siebert,
*Convergence of adaptive finite element methods*, SIAM Review 44(4), 631–658, 2002
(<https://doi.org/10.1137/S0036144502409093>); "MN" is K. Mekchay, R. H. Nochetto,
*Convergence of adaptive finite element methods for general second order linear elliptic PDEs*,
SIAM J. Numer. Anal. 43(5), 1803–1827, 2005 (<https://doi.org/10.1137/04060929X>).
Copies read (2026-10-17), in `../NavigatingIntrinsicTriangulations/cited/applications-5.4-5.5/`:
the authors' PDFs of VHM (19 pages) and GODF (8 pages); the preprint of MNS from Nochetto's
web page (27 pages, numbering of the preprint; its text layer is unreadable, pages 4–6 were
read as images); of MN only the first two pages of the authors' preprint, and in full
Mekchay's PhD thesis (University of Maryland 2005, `hdl:1903/3216`), whose Chapter 1 is the
paper; equation numbers below are the thesis's. Bibliographic data of the journal versions were
not checked against Crossref (no network lookup of the DOIs was made); the DOIs are those the
authors' pages give.

**Machine-checked** (`DdgProofs/Intrinsic/AppsAmrVector.lean`): `app-vector-face-holonomy`
(the transport angles around a face add up to `Σ (s − 1) θ` modulo full turns; the corner terms
at an interior vertex add up to its angle defect), `app-vector-fd-connection-laplacian` (one
edge's contribution to `u* L u` is `w |u_j − r u_i|²`; the entries are Hermitian),
`app-amr-map-distortion` (singular values from the conformal and anticonformal parts);
`app-vector-singularity-index` (`DdgProofs/Intrinsic/SingularityIndex.lean`, GODF Sec. 6.1.3,
Eq. 19, App. B): `ω_ij ∈ (−π, π]` (the range of `arg`/`atan2`) with
`u_j/|u_j| = e^{iω_ij} r_ij u_i/|u_i|` on every halfedge, except that at a tie (both halfedges of an
edge at `+π`) the halfedge with the higher index takes `−π`; on a closed valid oriented mesh with signpost frames, every vertex in some face,
every face index `(ω_ij + ω_jk + ω_ki + nΩ_ijk)/2π` is an integer and the indices add up to `nχ`
(discrete Poincaré–Hopf) unconditionally. (Evaluating `arg` on both halfedges, as the code did
before, gives `+π` twice at an edge with exactly opposite transported values, and the sum
`nχ + #such edges`.) Not
machine-checked: the entries below, and everything about rounding and the solvers.

### [formula: app-vector-transport] Transport of tangent vectors along an edge

- **Statement.** With `φ_ij` the direction of edge `ij` at `i` in the tangent frame of `i`, a
  vector is kept parallel from `i` to `j` by turning its angle by `ρ_ij = (φ_ji + π) − φ_ij`;
  `r_ij = e^{iρ_ij}`, and for `n`-vectors `r_ij^n`.
- **Reference.** VHM Sec. 5.2, "Discrete Parallel Transport", Eq. 10; GODF Sec. 6, Eq. 12 (the
  same with `θ_j(X_j, e_ij) − θ_i(X_i, e_ij)`).
- **Why not formalized.** It is the definition of the discrete connection (the Levi-Civita
  connection of the flat cone metric along a straight edge), not a derived formula. Its
  consistency around a face is machine-checked (`app-vector-face-holonomy`).
- **What depends on it.** Every matrix of `vector` (SPEC N40, N42).

### [formula: app-vector-fem-connection-laplacian] The finite element connection Laplacian of Knöppel et al.

- **Statement.** On a triangle with holonomy `Ω` (curvature constant per area) the piecewise
  linear basis sections `ψ_j = b_j ξ_j` (`ξ_j` parallel along rays from vertex `j`) have
  `⟨ψ_j, ψ_k⟩ = r̄_jk e^{iΩ b_i} b_j b_k` and the covariant derivatives of SPEC N43; their
  integrals are the local mass and Dirichlet matrices. In closed form
  `∫⟨ψ_j, ψ_k⟩ = r̄_jk |t| (6e^{iΩ} − 6 − 6iΩ + 3Ω² + iΩ³)/(3Ω⁴)`, `∫|ψ_i|² = |t|/6`, and the
  Dirichlet diagonal `(1/4|t|)[|p_jk|² + Ω² (3g₁₁ − 3g₁₂ + g₂₂)/90]`.
- **Reference.** GODF Sec. 6, "PL Finite Elements", Sec. 6.1.1 (Eq. 17, 18) and App. D.1–D.3.
- **Source check.** Read GODF Sec. 6 and App. D in the authors' PDF. The text layer of the
  closed forms of App. D.3 is garbled; the Rust code does not use them. It integrates the
  integrands derived in SPEC N43 exactly (monomial by monomial); the unit tests compare the
  mass entry with Eq. 17, the diagonal with the formula above, the flat limit with the cotan
  matrix, and all nine entries with a quadrature of the model (basis sections in one gauge,
  covariant derivatives by differences). The monomials are integrated with the area-coordinate
  rule `∫_T b_0^a b_1^b b_2^c = 2|T| a! b! c! / (a + b + c + 2)!` (e.g. O. C. Zienkiewicz,
  R. L. Taylor, *The Finite Element Method*, vol. 1); its case `a + b + c = 2` is the Galerkin
  mass matrix, machine-checked as `app-fem-galerkin-mass` and used by `amr` too.
- **Why not formalized.** The statement is about a gauge model of a line bundle that is not in
  Mathlib. Proved instead (`DdgProofs/TriangleIntegral/`, `[formula: fem-area-coordinate-rule]`):
  the area-coordinate rule as a Lebesgue integral over the triangle, and the slice formula that
  `monomial_integral` evaluates, `∫_T b₀^{e₀} b₁^{e₁} b₂^{e₂} e^{iΩ b_c} = 2|T| q! r!/(q + r + 1)! ·
  ∫₀¹ s^{e_c} (1 − s)^{q+r+1} e^{iΩs} ds` (`q = e_{c+1}`, `r = e_{c+2}`); and
  (`DdgProofs/TriangleIntegral/BetaExp.lean`, `[formula: app-vector-beta-exp]`) the two ways
  `beta_exp` evaluates the one-dimensional integral, exactly: the power series with its running
  product and a geometric bound on the tail after any term `K` with `|Ω| < K + 1`, and for `Ω ≠ 0`
  the recurrence `K_k = (e^{iΩ} − k K_{k−1})/(iΩ)` with the binomial expansion of `(1 − s)^m`.
  Not proved: the closed forms above, and that the integrands of SPEC N43 are those of the model.
- **What depends on it.** `ConnectionLaplacian::FiniteElement`: the smoothest direction fields
  and the energies of the README's Fig. 21 table.

### [formula: app-vector-maximum-principle] The vector maximum principle on an intrinsic Delaunay triangulation

- **Statement.** If every cotan weight `w_ij` at an interior vertex `i` is non-negative and
  `(L∇ u)_i = 0`, then `u_i = Σ_j w_ij r_ji u_j / Σ_j w_ij` is a positive combination of the
  transported neighbours and lies in their convex cone; on the intrinsic Delaunay
  triangulation every interior weight is non-negative.
- **Reference.** VHM Sec. 5.4; the paper Sec. 5.5 ("Flip-Free Vector Fields", Fig. 22).
- **Why not formalized.** Immediate from the definition (and from `it-delaunay-condition`
  for the weights); the margin of SPEC N45 is a measurement, not a theorem.
- **What depends on it.** `vector::check_maximum_principle` and the README's table.

### [formula: app-amr-residual-estimator] The residual a posteriori error estimator and the oscillation

- **Statement.** For `(c + Δ)u = f` and its linear finite element solution `u_H`, with the
  element residual `R_T = f − c u_H` (the Laplacian of a linear function vanishes inside an
  element) and the jump `J_S` of the normal derivative across side `S`,
  `η_T² = H_T² ‖R_T‖²_T + Σ_{S⊂∂T} H_S ‖J_S‖²_S`, `osc_T² = H_T² ‖R_T − R̄_T‖²_T` (`R̄_T`
  the mean), `H_T` the diameter of `T`; the energy error is bounded above by `C₁ η` and locally
  below by `C₂ η_T² − C₃ osc²` (upper and lower bounds).
- **Reference.** MN (thesis) Eqs. 1.3.3–1.3.8; MNS Sec. 2, Eqs. 2.5–2.8 (per side, with
  `‖H_S^{1/2} J_S‖²_S + ‖H f‖²_{ω_S}`); the paper's supplement, Sec. 2 (which prints
  `A_ijk²` in place of `H_T²`; SPEC N47).
- **Why not formalized.** Sobolev spaces and Clément interpolation are not in Mathlib in the
  needed form. Proved (`DdgProofs/TriangleIntegral/Fem.lean`, `[formula: app-amr-element-norms]`):
  the element norms as `estimate` evaluates them, `‖R‖²_T = (|T|/6)(Σ r_i² + Σ_{i<j} r_i r_j)` and
  `‖R − R̄‖²_T = ‖R‖²_T − |T| R̄²` for linear `R` with corner values `r`, `R̄` their mean (which is
  the mean of `R` over `T`); and the jump terms (`DdgProofs/TriangleIntegral/Flux.lean`,
  `[formula: app-amr-edge-flux]`): `side_flux` is `ℓ ∇u·ν` and the two fluxes of an interior side
  add up to `ℓ` times the jump of the normal derivative. The side term `ℓ_S ‖J_S‖²_S = (ℓ_S J_S)²`
  for a constant jump (a line integral) is not stated.
- **What depends on it.** `amr::estimate` and the marking; nothing of correctness of a solve.
  For a point source the estimator is not an upper bound of anything finite (SPEC N47).

### [formula: app-amr-bulk-marking] Bulk (Dörfler) marking

- **Statement.** Mark a set of elements carrying the fraction `θ` of the estimator
  (`Σ_marked η_T² ≥ θ² Σ η_T²`), and enlarge it until it carries the fraction `θ̂` of the
  oscillation; refining every marked element with an interior node reduces error plus
  oscillation by a fixed factor.
- **Reference.** MNS Sec. 3, Marking Strategy E (and O); MN Eqs. 1.3.9, 1.3.10 and Sec. 1.3.4
  (Interior Node Property); the supplement prints the criterion on sums of `η` (not squared)
  with `θ = 1/2` (SPEC N48).
- **Why not formalized.** The contraction is a theorem about Sobolev spaces; the marking itself
  is a greedy selection.
- **What depends on it.** `amr::mark`, `amr::adaptive_solve`. Not claimed: that intrinsic
  Delaunay refinement has the Interior Node Property (it inserts one circumcentre per marked
  face, not a node in every side), so the contraction is not claimed either (SPEC N48).

---

## Developability of Triangle Meshes

"Developability paper" is O. Stein, E. Grinspun, K. Crane, *Developability of Triangle Meshes*,
ACM TOG 37(4), Article 77, 2018 (<https://doi.org/10.1145/3197517.3201303>); section, equation
and proposition numbers are the paper's. The contract is `crates/ddg-developability/SPEC.md`
("SPEC", decisions `D1`–`D29`). Clean room: everything here is written from the paper's text and
the SPEC; the authors' code is never consulted. "L&O" is A. S. Lewis, M. L. Overton, *Nonsmooth
optimization via quasi-Newton methods*, Math. Program. A 141, 135–163, 2013 (numbering of the
preprint `cs.nyu.edu/~overton/papers/pdffiles/bfgs_inexactLS.pdf`, as the SPEC uses it); "N&W" is
J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer 2006.

**What the machine-checked part covers, and what it does not.** There is no certifiable minimizer:
the library runs a descent method on a nonconvex, nonsmooth energy. Machine-checked
(`DdgProofs/Developability/*.lean`, 30 `dev-` ids, statements in the README): the per-face
quantities and their derivatives along arbitrary differentiable motions of the corners, with the
corrections D2 (`1/(2A)` in Eq. 8) and D3 (squared length in Eq. 9); the adjoint
(back-propagation) form of the gradient; the three forms of Eq. 1 (D1); the gradients of the
active pieces of all six energies (for `E_λ` and the intrinsic width the derivative of `xᵀ A x`
at a fixed `x`, which is the derivative of the eigenvalue only through the cited
`dev-simple-eigenvalue-derivative`); the circumcenter formula and the site signs of D5; the
corrected log-map adjoint of Eq. 11 (D6) and the altitude derivative of B.5.1 (D7); one Jacobi
rotation; the acceptance conditions of the weak Wolfe line search as the oracle defines it; that
4-1 midpoint subdivision keeps normals, angles and vertex normals and creates hinges; the
small-angle classification of the remeshing; `E_λ` and the intrinsic energy of the regular cone.
**Not** machine-checked: the entries below; that the gradient of the active piece chosen at a tie
(SPEC §3.8, D10) is an element of the Clarke generalized gradient (L&O Sec. 5.1 practice); the
convergence of the Jacobi sweeps, of L-BFGS and of the flow; the closed forms of SPEC §10.3 for
`E_P`, `λ_max` and the intrinsic max energy (tests only); the topology of remeshing and cutting;
the Taylor series of the log-map coefficients near `φ = 0` and the overflow guard of the Jacobi
rotation (numerical devices); anything about rounding.

### [formula: dev-hinge-coplanar] An embedded vertex star is a hinge iff its normals are coplanar

- **Statement.** An embedded vertex star `St(i)` forms a hinge (Def. 3.1: its triangles split
  into two edge-connected regions of constant normal) if and only if all its triangle normals lie
  in a common plane through the origin.
- **Reference.** Developability paper Sec. 3.3, Theorem 3.3 and its proof (the cross product of
  two distinct normals of adjacent triangles is parallel to their common edge and to the normal of
  the plane, so a change of normal across more than two edges forces two edges from `i` to be
  parallel, which embeddedness excludes).
- **Why not formalized.** Embeddedness of a simplicial star (injectivity of the piecewise linear
  map of the star) is not set up in Mathlib or in `lean/proofs`; the argument is short but
  topological.
- **What depends on it.** The interpretation of all six energies as measures of developability
  (each vanishes exactly on coplanar normals, SPEC §1). No formula of the code: the energies are
  defined and differentiated without it, and the code does not check embeddedness (D28).
  Proved instead: on a star whose normals are orthogonal to one vector the objectives of Eq. 3,
  Eq. 6 and the intrinsic objectives vanish (`coplanar_energies_vanish`,
  `dev-midpoint-subdivision`).

### [formula: dev-hinge-antiparallel] A non-flat interior hinge has a pair of antiparallel edges

- **Statement.** If an interior vertex `i` is a non-flat hinge vertex, it is contained in a pair
  of antiparallel edges `ia, ib ∈ St(i)` (the two edges shared by the two flat regions, both
  parallel to `N₁ × N₂`).
- **Reference.** Developability paper App. A, Proposition A.1 and its proof.
- **Why not formalized.** Needs embeddedness (as `dev-hinge-coplanar`).
- **What depends on it.** Nothing in the code; it motivates `dev-hinge-ruled` and
  `dev-valence-three-flat`.

### [formula: dev-hinge-ruled] A discrete developable immersion without flat vertices is discrete ruled

- **Statement.** A discrete developable immersion `f` (every vertex star a hinge) with no flat
  vertices is discrete ruled (Def. 3.2): every vertex lies on a path of parallel edges with both
  endpoints on the boundary.
- **Reference.** Developability paper App. A, Proposition A.2 and its proof (by
  Proposition A.1, followed from vertex to vertex).
- **Why not formalized.** Needs `dev-hinge-antiparallel` and a global combinatorial argument on
  the mesh.
- **What depends on it.** Nothing in the code (it explains the rulings of the results and the
  edge weights of the cut, `dev-cut-edge-weight`).

### [formula: dev-valence-three-flat] A valence-3 hinge vertex is flat

- **Statement.** Any valence-3 hinge vertex is flat: if it were not, by Proposition A.1 it would
  have two antiparallel edges, which at valence 3 are edges of one triangle, so the triangle is
  degenerate and the star is not embedded.
- **Reference.** Developability paper App. A, Proposition A.3; Sec. 4.1.4 ("Triple points",
  Fig. 16), which uses it to leave valence-3 vertices out of the energy.
- **Why not formalized.** Needs `dev-hinge-antiparallel` (embeddedness).
- **What depends on it.** Only the motivation of the default `omit_valence_three = true`
  (SPEC §3.1, D8): a valence-3 interior vertex cannot become a non-flat hinge, so penalizing it
  would flatten it. Not a correctness property of any formula.

### [formula: dev-simple-eigenvalue-derivative] Derivative of a simple eigenvalue of a symmetric matrix

- **Statement.** If `λ` is a simple eigenvalue of the real symmetric matrix `A₀` with unit
  eigenvector `x₀`, there is a neighbourhood of `A₀` on which the eigenvalue `λ(A)` near `λ`
  is a smooth (real analytic) function of the symmetric matrix `A`, with differential
  `dλ = x₀ᵀ (dA) x₀`, i.e. `∇_A λ = x₀ x₀ᵀ`.
- **Reference.** J. R. Magnus, *On differentiating eigenvalues and eigenvectors*, Econometric
  Theory 1(2), 179–191, 1985, Theorem 1; T. Kato, *Perturbation Theory for Linear Operators*,
  Springer, 1995 reprint of the 2nd ed., Ch. II §5. Used by
  the developability paper App. B.3 ("the identity `∇_A λ = x xᵀ`").
- **Why not formalized.** Mathlib has the spectral theorem and the min-max (Rayleigh)
  characterization, but no perturbation theory of eigenvalues (differentiability of a simple
  eigenvalue as a function of the matrix).
- **What depends on it.** The gradient of `E_λ` (Eq. 10) and of the intrinsic width (SPEC §3.4,
  §3.6): proved are the Rayleigh characterization (`dev-covariance-matrix`), the supergradient
  inequality, which holds at every matrix (`dev-eigenvalue-supergradient`), and the derivative
  of `xᵀ A(t) x` at fixed `x` (`dev-covariance-gradient`, `dev-tangent-eigen-derivative`); this
  entry turns the last into the derivative of the eigenvalue. At a repeated smallest eigenvalue
  (every regular cone with `sin²α/2 < cos²α`) the energy is not differentiable and the code uses
  the gradient of one piece (SPEC §3.8); finite-difference tests are restricted to stars with an
  eigenvalue gap (SPEC §10.2).

### [formula: dev-max-covariance-voronoi] The minimum of `max_σ ⟨u, N_σ⟩²` is attained at a spherical Voronoi vertex

- **Statement.** For unit vectors `N_σ`, `φ(u) = max_σ ⟨u, N_σ⟩²` attains its minimum over the
  unit sphere at a point equidistant from three sites of `{±N_σ}` (a vertex of their spherical
  Voronoi diagram), or, when all normals are parallel or antiparallel, `φ ≡ 0` on a great circle;
  hence the minimum over the candidates of SPEC §3.5 (circumcenters of triples of sites, poles of
  pairs) is the minimum of `φ`.
- **Reference.** Developability paper App. B.4 (stated, with Fig. 28); the proof sketch of SPEC
  §3.5 (within a Voronoi cell the distance to the site has no interior maximum except on the
  great circle orthogonal to it, along a Voronoi edge `|⟨u, s⟩|` is monotone in the distance from
  the midpoint of the bisector), checked there against dense search on 200 random stars.
- **Why not formalized.** Spherical Voronoi diagrams and the largest empty circle problem are
  not in Mathlib; the sketch would need the cell decomposition of the sphere.
- **What depends on it.** That the code's `λ_max` is the minimum and not only an upper bound:
  every candidate is a unit vector, so the computed value is always `≥ λ_max` (exact
  arithmetic), and extra candidates are harmless. Proved: the circumcenter is equidistant from
  its three sites and is, up to sign, the only equidistant unit vector
  (`dev-spherical-circumcenter`), and the sensitivities at a generic minimizer
  (`dev-max-covariance-gradient`).

### [formula: dev-min-width-hull-edge] The minimum width of a convex polygon is attained orthogonally to an edge

- **Statement.** The width `max_σ ⟨P_σ, u⟩ − min_σ ⟨P_σ, u⟩` of a finite planar point set,
  minimized over unit directions `u`, attains its minimum at a direction orthogonal to an edge of
  the convex hull (rotating calipers).
- **Reference.** M. E. Houle, G. T. Toussaint, *Computing the width of a set*, IEEE Trans.
  Pattern Anal. Mach. Intell. 10(5), 761–765, 1988; G. T. Toussaint, *Solving geometric problems
  with the rotating calipers*, Proc. IEEE MELECON '83, Athens, 1983. Used by the developability
  paper App. B.5.1.
- **Why not formalized.** Needs convex hulls of planar point sets and the piecewise structure of
  the support function; not set up in `lean/proofs`.
- **What depends on it.** That the intrinsic max energy over the pair candidates
  `u_ab = n0 × (Ñ_b − Ñ_a)/|Ñ_b − Ñ_a|` (D7) is the minimum and not only an upper bound (every
  candidate is a unit tangent vector, so the computed value is always an upper bound). Proved:
  the value as `(max − min)²` (`dev-intrinsic-width`) and the derivatives
  (`dev-width-altitude-gradient`, `dev-min-width-gradient`).

### [formula: dev-lbfgs-two-loop] The L-BFGS two-loop recursion

- **Statement.** For stored pairs `(s_i, y_i)` with `s_iᵀ y_i > 0` and `H₀ = γ I`, `γ > 0`, the
  two-loop recursion returns `H_k g`, where `H_k` is the matrix of the limited-memory BFGS update
  of `H₀` by the stored pairs; `H_k` is symmetric positive definite, so `d = −H_k g` is a descent
  direction for `g ≠ 0`.
- **Reference.** N&W Algorithm 7.4 (two-loop recursion), Algorithm 7.5 (L-BFGS), Sec. 7.2,
  Eq. (7.19) (the scaling `γ = sᵀy/yᵀy`); positive definiteness of the BFGS update for
  `sᵀy > 0`: N&W Sec. 6.1.
- **Why not formalized.** Linear algebra of rank-two updates; the code uses the result only to
  expect a descent direction, which it checks at run time (`gᵀd < 0`, otherwise a reset to
  steepest descent, SPEC §4.2).
- **What depends on it.** The efficiency of the optimizer, not its safety: every accepted step
  satisfies the line-search conditions, which are machine-checked (`dev-lewis-overton-line-search`),
  and pairs with `sᵀy ≤ 0` are skipped (`dev-weak-wolfe-curvature` says when `sᵀy > 0` holds).

### [formula: dev-lewis-overton-termination] Existence of weak Wolfe steps and termination of the bracketing line search

- **Statement.** If `f` is locally Lipschitz and bounded below along the ray and `d` is a descent
  direction (`s = gᵀd < 0`), with `0 < c1 < c2 < 1`, there are steps satisfying the Armijo and
  weak Wolfe conditions, and the bracketing (doubling and bisection) search of Algorithm 2.6
  terminates for semi-algebraic `f` (or more generally under the paper's hypotheses on the
  one-dimensional function).
- **Reference.** L&O Theorem 2.5 (existence) and Theorem 2.7 (termination of Algorithm 2.6), as
  numbered in the SPEC (§4.3); the defaults `c1 = 0`, `c2 = 1/2` are those of L&O Sec. 5.1.
- **Source check.** The theorem numbers are taken from the SPEC; the preprint was not re-read for
  this entry.
- **Why not formalized.** Needs nonsmooth analysis (Clarke subdifferentials, semi-algebraic
  functions) not in Mathlib.
- **What depends on it.** Only that the line search usually succeeds within its limits; the code
  bounds the number of doublings and bisections and reports a weak or failed search
  (`LineSearchFailed`), and the machine-checked part covers what an accepted, capped or weak step
  satisfies (`dev-lewis-overton-line-search`).

### [formula: dev-armijo-backtracking] Finite termination of backtracking for a descent direction

- **Statement.** If `f` is continuously differentiable near `x`, `gᵀd < 0` and `0 < c1 < 1`, then
  `f(x + t d) ≤ f(x) + c1 t gᵀd` for all sufficiently small `t > 0`; hence backtracking
  `t ← ρ t` (`0 < ρ < 1`) stops after finitely many reductions.
- **Reference.** N&W Algorithm 3.1 (backtracking line search) and the discussion after it in
  Sec. 3.1; existence of steps satisfying the Wolfe conditions: N&W Lemma 3.1.
- **Why not formalized.** A one-line consequence of the definition of the derivative; the energy
  is only piecewise smooth, so the hypothesis holds only away from the switching set, and the
  code bounds the number of reductions anyway (`max_shrinks`).
- **What depends on it.** That the backtracking variant of gradient descent usually accepts a
  step; a failure ends the run with `LineSearchFailed` (SPEC §4.4).

### [formula: dev-flip-manifold] Flipping an interior edge keeps a manifold triangle mesh

- **Statement.** Replacing the two triangles `ijk`, `jil` of an interior edge `ij` of a manifold
  triangle mesh by `kli`, `lkj` gives a manifold triangle mesh with the same vertices, boundary
  and Euler characteristic, provided `k ≠ l` and `k, l` are not already joined by an edge.
- **Reference.** H. Edelsbrunner, *Geometry and Topology for Mesh Generation*, Cambridge
  University Press, 2001, Ch. 2 (edge flips), as cited by the SPEC; the same operation on
  surfaces is the flip of `it-delaunay-flip-termination` above (signpost paper Sec. 3.3).
- **Source check.** Not re-read for this entry; the statement is elementary combinatorics.
- **Proved (trimmed, 2026-11-05).** The topological part, for the flip on gluing data
  (`DdgProofs/Delaunay/Flip.lean`, `it-complex-flip`) and hence for an oriented mesh whose vertex
  stars are single fans (`DdgProofs/Delaunay/FromTriMesh.lean`, `it-complex-from-mesh`): the flip
  of an interior edge whose two triangles are distinct keeps the vertices (each with its class, also
  when `i, j, k, l` coincide), the number of edges and faces, the boundary halfedges with their end
  points, the Euler characteristic, the orientation, and the property that every vertex star is one
  fan (`Gluing.flipVertexEquiv`, `flip_vertexOf_surviving`, `flip_diamond_vertices`,
  `flip_eulerChar`, `flip_numBoundary`, `Faithful.flip`, `TriMesh.eulerChar_toGluing`).
- **Still cited.** That the flipped mesh is again *simplicial* in the sense of `TriMesh`
  (no face with a repeated vertex, no two edges with the same end points: this is where `k ≠ l`
  and "`kl` is not an edge" enter) and satisfies `TriMesh.Manifold` (every edge in at most two
  faces, every boundary vertex on exactly two boundary edges), i.e. that the remesher's re-listing
  of the two faces as `(a, q, p)`, `(q, b, p)` is the gluing-data flip up to the rotation of the
  faces. (That every boundary fan has exactly one boundary halfedge leaving and one entering its
  vertex is now proved in the gluing model, for any involution `twin`:
  `Gluing.boundaryCount_le_one`, `boundaryInCount_eq_boundaryCount`,
  `DdgProofs/Delaunay/BoundaryFan.lean`, `complex-boundary-fan`; what remains is the translation
  of the flipped gluing data back to a simplicial `TriMesh`.)
- **Why not formalized.** The remaining step is the simplicial bookkeeping between the two models;
  the Rust code rebuilds and validates the mesh after every remeshing pass (`TriMesh::new`, SPEC
  §5.3), so a violation is an error, never a silent corruption.
- **What depends on it.** That legal flips are never refused by the validation.

### [formula: dev-link-condition] An edge collapse keeps a 2-manifold iff the link condition holds

- **Statement.** Contracting an edge `pq` of a triangulated 2-manifold (with boundary) gives a
  homeomorphic triangulated 2-manifold if `link(p) ∩ link(q) = link(pq)` (and, at the boundary,
  the corresponding condition for boundary edges).
- **Reference.** T. K. Dey, H. Edelsbrunner, S. Guha, D. V. Nekhayev, *Topology preserving edge
  contraction*, Publ. Inst. Math. (Beograd) (N.S.) 66(80), 23–45, 1999 (the link condition for
  2-manifolds), as cited by the SPEC.
- **Source check.** Not re-read for this entry.
- **Why not formalized.** Simplicial topology of the mesh, not modelled in `lean/proofs`; as for
  flips, the result of every pass is rebuilt and validated (SPEC §5.3).
- **What depends on it.** That legal collapses are never refused by the validation and do not
  change the topology.

### [formula: dev-cut-edge-weight] Edge weights of the cut

- **Statement.** The cut through the high-energy vertices is a shortest network for the weights
  `w_e = c_len ℓ_e/ℓ̄ + c_align (1 − ½(|ê·x_a| + |ê·x_b|))` (`x_v` the eigenvector of the smallest
  eigenvalue of `A_v`; defaults `c_len = 1`, `c_align = 0`).
- **Reference.** Developability paper Sec. 4.4 ("a linear combination of length and the dot
  product between the edge direction and the smaller eigenvector"); SPEC §6 step 2, D22 (the
  combination and sign are a choice, open question 2 of the SPEC).
- **Why not formalized.** A design choice, not a derived formula; the weights are nonnegative by
  construction (`|ê·x| ≤ 1`).
- **What depends on it.** Where the seams go, not whether the cut is valid: every chart is a disk
  whatever the weights (guaranteed by `plan_cut`).

## Repulsive Curves

"RC paper" is C. Yu, H. Schumacher, K. Crane, *Repulsive Curves*, ACM TOG 40(2), Article 10, 2021
(<https://doi.org/10.1145/3439429>); section, equation, lemma and appendix numbers are the paper's.
The contract is `crates/ddg-repulsive-curves/SPEC.md` ("SPEC", decisions `D1`–`D30`, formula
table §12). "N&W" is J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer 2006.

**What the machine-checked part covers, and what it does not.** There is no certifiable minimizer:
the library runs a descent method on a nonconvex energy. Machine-checked
(`DdgProofs/RepulsiveCurves/*.lean`, 44 `rc-` ids, statements in the README): the edge quantities
and dual lengths; the kernel's identities, invariances and homogeneity, its gradient for `ν > 0`
and the continuity of that gradient at `ν = 0` for `α > 1` (D24); the discrete energy of Eq. 17–18,
its homogeneity `c^{2+α−β}`, invariance under translations, orthogonal maps and edge re-orientation,
the derivative of every pair term along any differentiable motion (SPEC §3.2), the assembled
vertex gradient, and the identities `Σ g = 0`, `Σ γ × g = 0`, `Σ γ·g = (2+α−β)Ê`; the Gram
matrices `B`, `B⁰` assembled by the increments of Sec. 5.2.2–5.2.3 and their quadratic forms
(Eq. 21), symmetry, positive semidefiniteness and kernels, the Eq. 31 identity as matrices with
`K = w + wᵀ`, Eq. 32 as the midpoint weights and the D2 kernel identity, the scaling of both
quadratures (`B, B⁰ ∝ c^{−(2σ+1)}`, D3), the H¹/H² comparison metrics; the saddle systems of
Eq. 22 and 23 (closest point, existence and uniqueness, descent, least-norm chord step), the
projector of D4 and Eqs. 30 and 33, the rows of §6.4, the scaling of the direction; the constraint
Jacobians of §5.1 with D6, D7 and the Duff et al. frame; the penalties of §5.3 (the surface
potential given the cited distance gradient); the length schedule and the subdivision mapping;
the conservative-advancement lemma of the CCD; the lumped quantities of the BVH, Eq. 26 for a
collapsed cluster, the B.1.2 increment, the rank-one BCT block; the coarsening rule and the
prolongation. **Not** machine-checked: the entries below; anything about the discrete energy as
an approximation of the smooth one; the derivative of the kernel at `ν = 0` (only the continuity
of the gradient there); the convergence of the flow, of the chord projection, of multigrid, of
L-BFGS and NCG; the error of Barnes–Hut and of the block cluster tree; that the BFS of the code
satisfies `CoarsenRule` (it is the rule written as a predicate); the data structures (BVH, BCT,
the multigrid hierarchy) and every tolerance; anything about rounding.

### [formula: rc-exponent-range] Blatt's lemma: the energy space of the tangent-point energy

- **Statement.** For `α > 1`, `β ∈ [α + 2, 2α + 1)` and `s = (β − 1)/α ∈ (1, 2)`, an embedded curve
  `γ ∈ C¹(S¹; R³)` has finite tangent-point energy `ℰ^α_β(γ)` if and only if, up to
  reparametrization, `γ ∈ W^{s,α}(S¹; R³)`.
- **Reference.** RC paper Lemma A.1 (App. A.2), citing S. Blatt, *The energy spaces of the tangent
  point energies*, J. Topol. Anal. 5(3), 261–270, 2013 (the main theorem), and S. Blatt,
  P. Reiter, *Regularity theory for tangent-point energies: the non-degenerate sub-critical case*,
  Adv. Calc. Var. 8(2), 93–116, 2015.
- **Source check.** Read in the RC paper; Blatt's article was not re-read for this entry.
- **Why not formalized.** Fractional Sobolev–Slobodeckij spaces on curves and the
  characterization of the energy space are not in Mathlib.
- **What depends on it.** Only the admissible exponent range: `Exponents::new` rejects `α ≤ 1`,
  `β < α + 2`, `β ≥ 2α + 1` (D24). The discrete formulas and their proofs hold for any real
  exponents.

### [formula: rc-differential-order] The differential of the tangent-point energy has order `2s`

- **Statement.** `dℰ^α_β` maps `W^{s,α}` to its dual `W^{−s,α'}`; in this sense it is a
  differential operator of order `2s = 2σ + 2`, which motivates an inner product of the same order
  (`B_σ` of Eq. 12, order `2σ + 2`; `B⁰_σ` of Eq. 13, lower order).
- **Reference.** RC paper App. A.3 ("Order of the Differential") and A.4.1; the dual space
  `W^{−s,q} = (W^{s,p})*`: E. Di Nezza, G. Palatucci, E. Valdinoci, *Hitchhiker's guide to the
  fractional Sobolev spaces*, Bull. Sci. Math. 136(5), 521–573, 2012, Remark 2.5 (as cited in
  App. A.1). A.3 writes `W^{s,p}` for `W^{s,α}` (SPEC D3).
- **Why not formalized.** A heuristic statement about smooth operators between fractional
  spaces; it motivates the choice of the metric and is not used by any formula.
- **What depends on it.** The efficiency of the Sobolev descent (the step sizes and iteration
  counts of Sec. 7), not its correctness: every accepted step is a descent step
  (`rc-descent-direction`) checked by the line search.

### [formula: rc-fractional-operator] The fractional Laplacian and `L_σ` are of order `2σ`

- **Statement.** For `σ ∈ (0, 1)` the fractional Laplacian `(−Δ)^σ` (Eq. 24) and the operator
  `L_σ` (Eq. 25, extrinsic distances `|γ(x) − γ(y)|` in place of intrinsic ones) are fractional
  operators of order `2σ`; `B_σ = 𝒟* L_σ 𝒟` has order `2σ + 2`.
- **Reference.** RC paper App. A.1.1 (Eq. 24), A.4 (Eq. 25) and A.4.1 (`⟨⟨B_σ u, v⟩⟩ =
  ⟨⟨L_σ 𝒟u, 𝒟v⟩⟩`); Di Nezza, Palatucci, Valdinoci 2012 (above), Sec. 3 (the fractional Laplacian
  as a singular integral).
- **Why not formalized.** Smooth theory of nonlocal operators, not in Mathlib.
- **What depends on it.** Nothing in the discrete formulas: the discrete metric is defined by
  Eq. 20–21 (and 31–32), whose algebra is machine-checked (`rc-high-order-gram`,
  `rc-low-order-gram`, `rc-kernel-matrix-form`, `rc-gram-psd-kernel`).

### [formula: rc-quadrature-consistency] Omitting adjacent edge pairs is insignificant under refinement

- **Statement.** The discrete energy of Eq. 17 omits pairs of edges with a common vertex; because
  the kernel is only weakly singular, these terms contribute nothing in the limit of regular
  refinement of a smooth curve, so Eq. 17 approximates `ℰ^α_β` of any sufficiently smooth curve
  through the vertices.
- **Reference.** RC paper Sec. 5.1 and Fig. 9; the energy of a polygon is infinite: P. Strzelecki,
  H. von der Mosel, *Geometric curvature energies: facts, trends, and open problems*, in *New
  Directions in Geometric and Applied Knot Theory*, De Gruyter 2017, Fig. 2.2 (as cited in
  Sec. 5.1).
- **Source check.** Read in the RC paper; Strzelecki–von der Mosel not re-read.
- **Why not formalized.** A consistency statement of a quadrature for a singular integral; the
  paper gives an argument, not a proof.
- **What depends on it.** That minimizers of the discrete energy resemble those of the smooth one.
  The flow minimizes the discrete energy as defined (`rc-energy`), whose formulas are checked.

### [formula: rc-chord-convergence] Local linear convergence of the chord method

- **Statement.** For `Φ` continuously differentiable near a regular point `γ*` of `{Φ = 0}`
  (`C = dΦ` of full row rank) and a starting point close enough to `γ*`, the iteration
  `γ̃ ← γ̃ + x`, `x` the `Ā`-least-norm solution of `C(γ) x = −Φ(γ̃)` with `C` and `Ā` frozen at a
  nearby `γ` (Eq. 23), converges linearly to a point of `{Φ = 0}`.
- **Reference.** J. M. Ortega, W. C. Rheinboldt, *Iterative Solution of Nonlinear Equations in
  Several Variables*, Academic Press 1970, Ch. 10 (Ostrowski's theorem on the convergence of
  fixed-point iterations, 10.1.3); C. T. Kelley, *Iterative Methods for Linear and Nonlinear
  Equations*, SIAM 1995, Ch. 5 (the chord method, Theorem 5.4.1), for the square case; the
  underdetermined case with the least-norm correction is a Gauss–Newton/chord iteration on the
  constraint manifold.
- **Source check.** The references are the SPEC's (§12); the theorem numbers were not re-read.
- **Why not formalized.** Needs a contraction argument with the inverse function theorem in
  several variables; the code needs only what one iteration does, which is machine-checked
  (`rc-chord-projection`).
- **What depends on it.** That the projection usually succeeds within 3 iterations. The code checks
  every row against its tolerance and rejects the trial step otherwise (SPEC §6.2), so a failure is
  reported, never accepted.

### [formula: rc-distance-gradient] The distance to a closed convex set is `C¹` outside it

- **Statement.** For a nonempty closed convex set `K ⊂ R³` (a triangle) with metric projection
  `π`, `½ δ²` with `δ(x) = |x − π(x)|` is differentiable with gradient `x − π(x)`; hence `δ` is
  differentiable at every `x ∉ K` with `∇δ = (x − π(x))/δ(x)`.
- **Reference.** H. H. Bauschke, P. L. Combettes, *Convex Analysis and Monotone Operator Theory in
  Hilbert Spaces*, 2nd ed., Springer 2017, Ch. 12 (Moreau envelopes; the Moreau envelope of the
  indicator of `K` is `½δ²`, its gradient is `Id − P_K`); the closest point on a triangle:
  C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann 2004, §5.1.5.
- **Source check.** As cited by the SPEC; the proposition numbers were not re-read.
- **Why not formalized.** Mathlib has the orthogonal projection onto closed subspaces and convex
  sets but not the differentiability of the distance function in this form.
- **What depends on it.** The gradient of the surface potential (D19); its chain rule given `∇δ`
  is machine-checked (`rc-penalty-surface`), and the Rust test compares it with finite differences.

### [formula: rc-near-cos-series] The two series of `S(θ) = ∫_0^θ cos^m φ dφ` (WP16, SPEC §3.4)

- **Statement.** For `m > 0` and `|θ| ≤ π/4`, with `s = sin θ`, `a = (m − 1)/2`:
  `S(θ) = Σ_{k≥0} C(a, k)(−1)^k s^{2k+1}/(2k + 1)`; for `ψ = π/2 − |θ| ≤ π/4`, `c = cos θ = sin ψ`:
  `∫_0^ψ sin^m = c^{m+1} Σ_{k≥0} C(−½, k)(−1)^k c^{2k}/(m + 1 + 2k)`, and
  `S(θ) = ±(S(π/2) − ∫_0^ψ sin^m)`. Both follow from the substitution `t = sin φ`
  (`∫_0^{sin θ} (1 − t²)^{(m−1)/2} dt`, resp. `∫_0^{sin ψ} t^m (1 − t²)^{−1/2} dt`) and termwise
  integration of the binomial series `(1 − x)^a = Σ C(a, k)(−x)^k`, which converges absolutely and
  uniformly for `|x| ≤ ½` (here `x = s²` or `c²` `≤ ½`). Equivalently
  `S(θ) = ½ B_{sin²θ}(½, (m+1)/2)` (incomplete beta function) with its hypergeometric series.
  The crate sums 60 terms; with ratio at most ½ the remainder is below `2⁻⁶⁰` times the leading
  term (`|C(a, k)|` grows at most polynomially in `k`).
- **Reference.** NIST *Digital Library of Mathematical Functions* (<https://dlmf.nist.gov>): §8.17
  (incomplete beta function `B_x(a, b)`, eq. 8.17.1; hypergeometric representation
  `B_x(a, b) = (x^a/a) F(a, 1 − b; a + 1; x)`, eq. 8.17.7; with `a = ½`, `b = (m + 1)/2` this is the
  second form of the statement) and §15.2 (the Gauss series, eq. 15.2.1, defined by the series on
  `|z| < 1`); the binomial series with general exponent: M. Abramowitz, I. A. Stegun, *Handbook of
  Mathematical Functions*, NBS 1964, §3.6 (binomial series, pp. 14–15 of the Dover printing; the
  equation number 3.6.8 of the first version of this entry was not verified and is dropped).
- **Source check.** Checked 2026-10-27 against <https://dlmf.nist.gov/8.17> and
  <https://dlmf.nist.gov/15.2>: eqs. 8.17.1, 8.17.7 and 15.2.1 read as above (8.17.8 is the
  different form with the factor `(1 − x)^b`). The derivation in the Statement (substitution
  `t = sin φ`, termwise integration, `S = ½ B_{sin²θ}(½, (m+1)/2)`) was rechecked by hand and is
  consistent with 8.17.7. A&S equation number not verified (no readable scan reached). The series are checked against
  adaptive Gauss–Legendre quadrature of `cos^m` to `2e-15` absolute for three exponents
  (`tests/near_pair.rs`, `series_and_inner_integral_match_brute_force`) and against SciPy to
  `4e-16` (WP16 prototype).
- **Why not formalized.** Termwise integration of a power series and the binomial series with real
  exponents are in Mathlib only in pieces; the identity is classical and the numbers are checked.
  What is machine-checked is that any `S` with `S' = cos^m` gives the inner integral
  (`rc-near-inner`, `near_inner_hasDerivAt`, `near_inner_integral`), and that `∫_0^θ cos^m` is one.
- **What depends on it.** The value of the inner integral `F(y)` (SPEC §3.4); its differential uses
  `S' = cos^m`, exact up to the truncation of the series.

### [formula: rc-near-closest] Closest points of two segments and their derivatives (WP16, SPEC §3.4)

- **Statement.** The closest points `P(s*)`, `Q(t*)` of two non-intersecting segments are found by
  minimizing the convex quadratic `|P(s) − Q(t)|²` over `[0, 1]²` (clamping, the parallel case by a
  choice). The distance `d = min |P(s) − Q(t)|` has the derivative
  `∇d = n̂ ⊗ ((1 − s*), s*, −(1 − t*), −t*)` (`n̂` the unit vector between the closest points)
  wherever the minimizer is unique (Danskin's theorem: the constraint set `[0, 1]²` does not depend
  on the positions); a free parameter (`0 < s* < 1`, resp. `t*`) satisfies `A·D = 0` (`B·D = 0`),
  and the implicit function theorem gives `∇s*`, `∇t*` from the `2 × 2` (or `1 × 1`) system of
  those conditions; clamped parameters have zero derivative. `d` is `C¹` but not `C²` where the
  set of active bounds changes, and not differentiable for parallel overlapping segments (the
  minimizer is not unique).
- **Reference.** C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann 2004, §5.1.9
  "Closest Points of Two Line Segments" (about p. 148; page not verified);
  J. M. Danskin, "The theory of max-min, with applications", *SIAM J. Appl. Math.* 14 (4) (1966)
  641–664, and the monograph *The Theory of Max-Min and its Application to Weapons Allocation
  Problems*, Springer 1967 (Ökonometrie und Unternehmensforschung 5); J. F. Bonnans, A. Shapiro,
  *Perturbation Analysis of Optimization Problems*, Springer 2000, Chapter 4 (first-order
  differentiability of optimal value functions; section and theorem numbers not verified).
- **Source check.** Checked 2026-10-27: the Ericson section number and title 5.1.9 (and 5.1.5
  "Closest Point on Triangle to Point") against the Routledge table of contents
  (<https://www.routledge.com/Real-Time-Collision-Detection/Ericson/p/book/9781558607323>); the
  Danskin 1967 monograph (Springer, Berlin, series volume 5, viii + 126 pp.) against library
  records (<https://www.springer.com/gp/book/9783642460944>). The 1966 SIAM paper is added from
  memory of the standard citation for Danskin's theorem and was not re-read. Bonnans–Shapiro title,
  publisher and year confirmed by catalogue records; the section number `§4.1` of the first version
  was not verified and is replaced by the chapter. Danskin's theorem is used in its standard form
  (the minimum over a compact set, here `[0, 1]²`, of a function `C¹` in the parameter; its
  directional derivative is the minimum of the partials over the minimizers, hence the gradient at
  a unique minimizer). The Rust differential through
  these derivatives is checked by central differences (`tests/near_pair.rs`).
- **Why not formalized.** The case analysis of the algorithm and the envelope theorem on a
  rectangle; the consequence (the differential of the blended energy) is tested.
- **What depends on it.** The blend weight `w(d/ℓ_IJ)` and the breakpoints of the outer
  quadrature (SPEC §3.4).

### [formula: rc-near-quadrature-error] Accuracy of the graded outer quadrature (WP16, SPEC §3.4)

- **Statement.** For an integrand analytic in the Bernstein ellipse `E_ρ` of an interval and bounded
  there by `M`, the `N`-point Gauss–Legendre error (`N ≥ 2`) is at most
  `64 M / (15 (ρ² − 1) ρ^{2(N−1)})` times the half-length (equivalently `ρ^{−2n}` for `n + 1`
  points); after the substitution `t = c ± δ sinh w` a near singularity at distance `δ`
  from `c` lies at imaginary distance about `π/2` from the real `w` axis for every ratio `δ/L`,
  so a fixed number of panels in `w` gives a uniform relative accuracy as `δ/L → 0`.
- **Reference.** L. N. Trefethen, *Approximation Theory and Approximation Practice*, SIAM 2013
  (Extended Edition 2019), Chapter 19, Thm. 19.3 (Gauss and Clenshaw–Curtis error for functions
  analytic in a Bernstein ellipse; the book's own page not seen); P. Rabinowitz, "Rough and ready
  error estimates in Gaussian integration of analytic functions", *Comm. ACM* 12 (5) (1969)
  268–270, ineq. (18) (the origin of the `64/15` constant); P. R. Johnston, D. Elliott, "A sinh
  transformation for evaluating nearly singular boundary element integrals", *Int. J. Numer.
  Meth. Engng.* 62 (4) (2005) 564–578, doi:10.1002/nme.1208 (the repository record of the Griffith
  copy gives the title with "two-dimensional" inserted; the wording of the title is as in the
  University of Tasmania eCite record).
- **Source check.** Checked 2026-10-27, secondarily, since the book text could not be read:
  Thm. 19.3 is cited for Gauss–Legendre quadrature of functions analytic and bounded by `M` in
  `E_ρ` in L. N. Trefethen, "Exactness of quadrature formulas" (SIAM Rev., 2022; volume and pages not verified; the passage is near p. 148,
  <https://people.maths.ox.ac.uk/trefethen/exactness.pdf>), and stated with the constant `64/15`
  for `N` points, `N ≥ 2`, and `ρ^{−2(N−1)}` in arXiv:2411.19086, Thm. 6.9
  (<https://arxiv.org/pdf/2411.19086>), which cites Rabinowitz ineq. (18) and Thm. 19.3 of the
  Extended Edition. Another paper (arXiv:2211.00696, Thm. 3.1) quotes Thm. 19.3 with `144/35`
  for `n + 1` points, `n ≥ 2`; so the constant (`64/15` vs `144/35`) and the exact hypotheses of
  the book's theorem are not verified. The first version of this entry wrote `ρ^{2n}` for `n`
  points, an off-by-one in the point count, now corrected. The Johnston–Elliott bibliographic data
  (journal, volume, issue, pages, year, doi) match the eCite record
  (<https://ecite.utas.edu.au/36544>, via search summary; the page itself was not readable); the
  claim "imaginary distance about `π/2`" was rechecked by hand (`sinh w = ±i` at `w = ±iπ/2`).
  The crate does not rely on the bound: the accuracy is
  measured (`tests/near_pair.rs`: at most `2.5e-6` relative against an adaptive reference on seven
  crafted pairs at distances `0.5 … 1e-6`, both exponent pairs; SPEC §17o).
- **Why not formalized.** An error bound of a numerical rule; only the substitution identity is
  machine-checked (`rc-near-outer`).
- **What depends on it.** How close the near-pair energy is to the pair integral of the polygon.

### [formula: rc-armijo] Finite termination of backtracking for a descent direction

- **Statement.** If `F` is continuously differentiable near `γ` and `dF·d < 0`, then
  `F(γ + τd) ≤ F(γ) + c τ dF·d` for all sufficiently small `τ > 0` (`0 < c < 1`), so halving
  `τ` stops after finitely many steps.
- **Reference.** N&W Sec. 3.1 and Algorithm 3.1 (backtracking line search); S. Boyd,
  L. Vandenberghe, *Convex Optimization*, Cambridge University Press 2004, Algorithm 9.2 (the
  paper's choice, Sec. 5.4). The same fact for the developability flow:
  `dev-armijo-backtracking`.
- **Why not formalized.** A consequence of the definition of the derivative; with the projection
  onto the constraint set and the CCD bound the accepted step is not the plain Armijo step, and
  the code bounds the number of halvings anyway.
- **What depends on it.** That the line search usually accepts a step; a failure ends the run with
  `Stuck` (SPEC §7.1–§7.3). That the direction is a descent direction is machine-checked
  (`rc-descent-direction`).

### [formula: rc-accd] The additive CCD algorithm

- **Statement.** Conservative advancement with the additive (Lipschitz) distance bound and a gap
  `s d₀` computes a conservative time of impact for linearly moving point–triangle and
  edge–edge pairs, in finitely many advancements for a positive gap.
- **Reference.** M. Li, D. M. Kaufman, C. Jiang, *Codimensional Incremental Potential Contact*,
  ACM TOG 40(4), Article 170, 2021, Sec. 5 (Additive CCD), as cited by the SPEC (§7.4).
- **Source check.** As cited by the SPEC; not re-read for this entry.
- **Why not formalized.** The algorithm's termination and efficiency; the safety of one
  advancement (no contact before the advanced time) is machine-checked
  (`rc-ccd-advancement`), and the code caps the number of advancements at 100, returning the
  current (still conservative) time.
- **What depends on it.** The efficiency of the collision-guarded line search (D11), not its
  safety.

### [formula: rc-isotopy] A collision-free path of embeddings keeps the isotopy class

- **Statement.** If every accepted step moves the curve network along a path (piecewise linear in
  time: the step and each chord iteration) along which no two edges without a common vertex
  intersect, and the edges keep positive length, then the network stays embedded along the path,
  and the initial and final networks are isotopic (and, being polygonal, ambient isotopic).
- **Reference.** RC paper Sec. 5.4 (the first time of collision as the start of the line search
  "guarantees that the curve remains in the same isotopy class"); for polygonal knots,
  G. Burde, H. Zieschang, *Knots*, 2nd ed., de Gruyter 2003, Ch. 1 (tame knots, isotopy).
- **Source check.** The textbook reference is generic and was not re-read for this entry.
- **Why not formalized.** Knot theory (isotopy of embeddings) is not in Mathlib. Note that
  adjacent edges are not checked by the CCD; a fold of two edges at a common vertex onto each other
  would need the angle between them to reach zero, which the energy (`k → ∞` when edges become
  close) and the length constraints resist but no theorem here excludes.
- **What depends on it.** The claim that the collision-guarded flow does not change the knot type
  (SPEC §7.2); the normalized Armijo search (D12) has no such guarantee.

### [formula: rc-bh-admissibility] The error of the Barnes–Hut approximation (Eq. 27)

- **Statement.** If `r_x^ℬ/|x_I − x̄_ℬ| ≲ ε` and `r_T^ℬ ≲ ε`, the lumped term of Eq. 26 approximates
  the sum over the node `ℬ` with an error of order `ε` relative to it (a Taylor expansion of the
  kernel about the lumped tangent-point).
- **Reference.** RC paper App. B.1.2 ("Admissibility", Eq. 27); J. Barnes, P. Hut, *A hierarchical
  O(N log N) force-calculation algorithm*, Nature 324, 446–449, 1986.
- **Why not formalized.** The paper gives a Taylor-series argument, not a quantitative bound; the
  exactness for a collapsed cluster and the derivative of the lumped term are machine-checked
  (`rc-bvh-lumped`, `rc-bh-differential`).
- **What depends on it.** The accuracy of the accelerated energy and differential; the SPEC
  measures it (§13.4: convergence to the exact values as `θ → 0`, equality at `θ = 0`). The
  endpoint-box guard that makes the exclusion of adjacent pairs exact is machine-checked
  (`guard_disjoint`, under `rc-bvh-lumped`).

### [formula: rc-bct-admissibility] Low-rank approximation of admissible blocks (Eq. 28)

- **Statement.** For an asymptotically smooth kernel, a block `(𝒜, ℬ)` with
  `max(r_x^𝒜, r_x^ℬ) ≲ ε |x̄_𝒜 − x̄_ℬ|` and `max(r_T^𝒜, r_T^ℬ) ≲ ε` is approximated by the rank-one
  matrix `ℓ[𝒜] κ(𝒜, ℬ) ℓ[ℬ]ᵀ` with an error of order `ε` (0th-order multipole expansion).
- **Reference.** RC paper App. B.2 (Eq. 28); W. Hackbusch, *Hierarchical Matrices: Algorithms and
  Analysis*, Springer Series in Computational Mathematics 49, 2015 (admissibility conditions and
  low-rank approximation of asymptotically smooth kernels); L. Greengard, V. Rokhlin, *A new version
  of the fast multipole method for the Laplace equation in three dimensions*, Acta Numerica 6,
  229–269, 1997.
- **Source check.** Chapter numbers of Hackbusch's book were not re-read for this entry.
- **Why not formalized.** Approximation theory of kernel functions, not in Mathlib; the algebra of
  the rank-one block and its exactness for collapsed clusters are machine-checked
  (`rc-bct-lowrank`).
- **What depends on it.** The quality of the multigrid preconditioner, not the correctness of the
  result: the residual of every solve is measured against the hierarchical product and reported
  (SPEC §8.4).

### [formula: rc-multigrid] Convergence of the V-cycle and of flexible CG for the projected system

- **Statement.** The multigrid V-cycle with conjugate-gradient smoothing on `PᵀĀPy = Pᵀa`
  (Eq. 30), with the coarsening and prolongation of B.3.2, reduces the residual by a factor bounded
  away from 1 independently of the resolution; flexible CG with the V-cycle as a variable
  preconditioner converges.
- **Reference.** RC paper App. B.3 (empirical: "6 or fewer V-cycles" to `10⁻³`); D. Braess,
  R. Sarazin, *An efficient smoother for the Stokes problem*, Appl. Numer. Math. 23(1), 3–19, 1997
  (multigrid for saddle-point problems via the projected system); M. Ainsworth, C. Glusa, *Aspects
  of an adaptive finite element method for the fractional Laplacian*, Comput. Methods Appl. Mech.
  Eng. 327, 4–35, 2017 (multigrid for fractional operators); Y. Notay, *Flexible conjugate
  gradients*, SIAM J. Sci. Comput. 22(4), 1444–1460, 2000.
- **Why not formalized.** There is no convergence proof for multigrid on `L_σ`-type operators with
  this coarsening (SPEC §12: "no proof for `L_σ`-type operators"); the algebra it relies on is
  machine-checked (`rc-projector`, `rc-projected-system`, `rc-projection-transform`,
  `rc-prolongation`, `rc-coarsen`).
- **What depends on it.** Speed. The solver stops on its measured residual or after 20 cycles and
  reports the residual history; stagnation switches to flexible CG (SPEC §8.4).

### [formula: rc-lbfgs-two-loop] The L-BFGS two-loop recursion

- **Statement.** The two-loop recursion with stored pairs `(s_i, y_i)`, `s_iᵀy_i > 0`, and an initial
  operator `H₀` (here the solve of the §6.4 system, "dynamic" metric) returns `H_k r`, where `H_k`
  is the L-BFGS update of `H₀`; `H_k` is symmetric positive definite when `H₀` is, so `−H_k r` is a
  descent direction.
- **Reference.** D. C. Liu, J. Nocedal, *On the limited memory BFGS method for large scale
  optimization*, Math. Program. 45, 503–528, 1989; N&W Algorithm 7.4 and Sec. 7.2. The same entry
  for the developability flow: `dev-lbfgs-two-loop`.
- **Why not formalized.** As for `dev-lbfgs-two-loop` (rank-two updates). With a positive definite
  `H₀` on `ker M_c` (machine-checked: `massRows_posDefOn`, `rc-mean-free-direction`).
- **What depends on it.** The efficiency of the evaluation method of Sec. 7 (soft problems only);
  pairs with `sᵀy ≤ 1e-12 |s||y|` are skipped (D22).

### [formula: rc-ncg-polak-ribiere] Polak–Ribière+ conjugate gradients with restarts

- **Statement.** With `β = max(0, ⟨r_k, z_k − z_{k−1}⟩/⟨r_{k−1}, z_{k−1}⟩)` (PR+) and the restart
  `d_k = −z_k` whenever `dF·d_k ≥ 0`, every direction used is a descent direction (by the restart
  rule) and the method converges globally under the usual line-search conditions.
- **Reference.** E. Polak, G. Ribière, *Note sur la convergence de méthodes de directions
  conjuguées*, Rev. Française Informat. Recherche Opérationnelle 3(16), 35–43, 1969; N&W Sec. 5.2
  (the PR+ choice, Eq. 5.45, and restarts).
- **Source check.** Equation number as in N&W 2nd ed., not re-read for this entry.
- **Why not formalized.** Global convergence of nonlinear CG is outside the scope of this
  package; that the restart yields a descent direction is immediate from the rule, and the
  preconditioned direction `z` itself is a descent direction (`rc-descent-direction`).
- **What depends on it.** The efficiency of the NCG evaluation method (soft problems only, D22).

## Integer Coordinates

"IC paper" is M. Gillespie, N. Sharp, K. Crane, *Integer Coordinates for Intrinsic Geometry
Processing*, ACM TOG 40(6), Article 252, 2021 (<https://doi.org/10.1145/3478513.3480522>; the
13-page paper with the supplement, Appendices B–D, in `../IntegerCoordinates/integer_coordinates.pdf`;
equation, algorithm and appendix numbers are the paper's). The contract is
`crates/ddg-integer-coords/SPEC.md` ("SPEC", decisions `I1`–`I26`, theorem table §9.2).

**What the machine-checked part covers, and what it does not.** Machine-checked
(`DdgProofs/IntegerCoords/*.lean`; the design of the combinatorial model is
`DdgProofs/IntegerCoords/DESIGN.md`): the arithmetic of Eq. 1, 2, 6 and of crossings; the **face
model** (one triangle of `T¹` with the pieces of the curves of `T⁰` in it, disjointness as
non-interleaving on the boundary cycle): every valid model is the canonical pairing of Fig. 11 /
Alg. 9 for corner counts given by Eq. 1–2, the side counts determine it, a model exists exactly for
the counts the oracle calls valid, at most one corner emanates, Alg. 1's step follows the arc; the
face split of I15 (the oracle's counts for a region are the numbers of arcs separating it from the
corners, and the new faces are valid for every region); Eq. 3 applied twice returns `n_ij` for every
pair of valid faces (via a quadrilateral model whose arcs are counted by kind); Thm. D.1 in the form
of its angle argument (and the failure of its strict form); chords of a convex quadrilateral and of
a triangle cross iff their ends interleave, cevians cross; App. A's recovery of positions for any
solution of Eq. 15; Eq. 5; Eq. 14 as the minimizer of Eq. 13 and the finite part of the mass
identity. Tier C (`DdgProofs/IntegerCoords/Complex/`, work package WP9): the **complex model**
(faces glued by a side pairing, Δ-complexes allowed, one face model per face, crossings glued by
`p ↦ n − 1 − p`, curves as the orbits of Alg. 1's step), `ic-trace-terminates`, and the labels
at the input vertices with `ic-identify-edge` (Alg. 3 computes the label of the arc; the walk back
past inserted vertices, Alg. 3 lines 1–5, is not modelled), the flip of a complex model with the
correctness of Eq. 3 after gluing two face models (`ic-flip-normal`, C part: same curves through the
quadrilateral), `ic-roundabout-update`, the combinatorial part of `ic-extract-intrinsic-edge` (each
crossing on exactly one curve from a vertex, found by walking back), and the merge of two halves (`ic-merge-edge-normal`: the
merged faces are valid and the Rust formula counts the merged edge's crossings, given the curve ends
at the inserted vertex as the realization has them); tier B: the split of a crossed edge
(`ic-split-crossed-edge`) and the split exactly at a crossing (`ic-split-at-crossing`: the new edge's
count as the arcs crossing it, the two new faces valid, the merge back); tier A: the strip intersection of Alg. 2 (`ic-strip-intersection`); and the
input-aligned frames (`ic-input-aligned-frame`, `DdgProofs/IntegerCoords/Frames/`: for the straight
segment of the strip, the Rust direction is the signposts' frame angle, at input and inserted
vertices; that the geodesic edge is that segment is `it-straightest-geodesic` with
`ic-normal-isotopy`). All statements
about integer formulas are about the oracle's own definitions
(`lean/oracle/DdgOracle/IntegerCoords/Formulas.lean`). **Not** machine-checked: the entries below;
the data structure of the Rust crate and that
it applies the formulas to the right elements; anything about rounding (which region floating point
picks, I15, is outside the proofs by design: every region is proved to give a valid state).

Sources were checked on 2026-10-31. Bibliographic data from the Crossref records of the DOIs where
one exists (Haken, Schaefer et al., Erickson–Nayyeri, Xia, Sharp–Crane 2020a and 2020b); Kneser 1929
has no DOI and is cited as the IC paper cites it. "Read" means the text was read; the IC paper was
read in full for these entries; the works it cites were not re-read unless the entry says so. No
formalization of any statement below was found in Lean/Mathlib, Isabelle (AFP), Coq, HOL Light or
Mizar (survey of 2026-10-31, `runs_tmp/ic-formalization-survey.md`; items of 2026 from that survey
are not cited here).

### [formula: ic-normal-isotopy] Normal coordinates determine the curves up to normal isotopy

- **Statement.** On a triangulated surface, a family of disjoint simple arcs between vertices that
  meet every triangle in normal arcs is determined, up to an isotopy that keeps every arc normal, by
  its normal coordinates (the number of crossings per edge); conversely, every vector of
  nonnegative integers satisfying the local conditions of every triangle (Eq. 1 integral and
  nonnegative) is the vector of such a family. Hence the integers `n` (with the roundabouts, which
  name the curves) encode `T⁰` over `T¹` (IC paper Sec. 2.3; SPEC §1, §2.4 item 4).
- **Reference.** H. Kneser, *Geschlossene Flächen in dreidimensionalen Mannigfaltigkeiten*,
  Jahresber. Dtsch. Math.-Ver. 38, 1929; W. Haken, *Theorie der Normalflächen: Ein
  Isotopiekriterium für den Kreisknoten*, Acta Math. 105, 245–375, 1961,
  <https://doi.org/10.1007/BF02559591> (normal surfaces, of which normal curves are the
  two-dimensional case); for curves on surfaces and their algorithms: M. Schaefer, E. Sedgwick,
  D. Štefankovič, *Algorithms for Normal Curves and Surfaces*, COCOON 2002, LNCS 2387, 370–380,
  <https://doi.org/10.1007/3-540-45655-4_40>; J. Erickson, A. Nayyeri, *Tracing Compressed Curves in
  Triangulated Surfaces*, Discrete & Computational Geometry 49(4), 823–863, 2013,
  <https://doi.org/10.1007/s00454-013-9515-z>, Sec. 2.3 (corner coordinates dual to normal
  coordinates, cited by the IC paper after Eq. 2). IC paper Sec. 2.3, 3.1.
- **Source check.** Read: IC paper Sec. 2.3 and its citations (Kneser 1929; Haken 1961; Schaefer et
  al. 2002; Erickson–Nayyeri 2013, Sec. 2.3). Bibliographic data of Haken, Schaefer et al.,
  Erickson–Nayyeri from Crossref. Not read: the four works themselves.
- **Why not formalized.** It is a statement about isotopy classes of curves on a topological
  surface; neither Δ-complexes, surfaces as gluings, nor isotopy of curves are in Mathlib, and no
  Jordan–Schoenflies theorem is available. What is proved is its finite combinatorial side for one
  triangle: the face counts determine the face model and exist exactly for valid counts
  (`FaceModel.Valid.partner_eq`, `FaceModel.exists_valid_iff`, `[formula: ic-face-decoding]`).
- **What depends on it.** That `n` and `r` are a complete encoding of the correspondence; the
  library's exactness claim (SPEC §1) is relative to it; so is the strip of input faces from which
  the input-aligned frames are read (`ic-input-aligned-frame`, proved for that strip). Tested by the certificate (SPEC §2.4 item 4)
  and the brute-force recount (SPEC §8 layer 3).

### [formula: ic-refinement-termination] Intrinsic Delaunay refinement terminates

- **Statement.** On a closed intrinsic triangulation whose vertex angle sums are at least `60°`,
  Alg. 10 (Chew's second algorithm on the intrinsic triangulation, circumcentres inserted by tracing,
  `T¹` kept Delaunay) terminates with all corner angles at least `30°` (IC paper Thm. C.1, with
  Lemmas C.2–C.5). With boundary (split encroached boundary edges, remove inserted vertices in the
  Dijkstra ball) the paper proves nothing and reports success on its dataset (Sec. 6.1, Table 1).
  Removing more inserted interior vertices than the ball requires does not affect termination
  (App. C.1, citing Shewchuk 1997, Sec. 3.4.2).
- **Reference.** IC paper App. C.2, Thm. C.1, Lemmas C.2–C.5 (C.2 rests on Indermitte et al. 2001,
  Proposition 1, C.3 and C.4 on Bobenko–Springborn 2007, Definition 3); the planar argument of
  J. R. Shewchuk, *Delaunay Refinement Mesh Generation*, PhD thesis, CMU-CS-97-137, 1997,
  Sec. 3.2.1 and Lemma 12; L. P. Chew, *Guaranteed-quality mesh generation for curved surfaces*,
  SoCG 1993 (see `it-chew-second-algorithm` above for both).
- **Source check.** Read: IC paper App. C (Thm. C.1 and its proof, Lemmas C.2–C.5, C.1). The planar
  sources as in the entries `it-chew-second-algorithm`, `it-ruppert-encroachment` above. Note: the
  proof of Lemma C.5 (sliding a geodesic loop to a vertex without changing its length) is informal
  in the paper.
- **Why not formalized.** It needs geodesics, Delaunay tessellations and packing arguments on
  piecewise flat surfaces (none in Mathlib); the boundary case is open in the paper itself.
- **What depends on it.** That refinement stops by itself; the library bounds it by an insertion
  limit and reports reaching it (SPEC §5, I23).

### [formula: ic-dijkstra-ball] The geodesic ball lies inside the Dijkstra ball of twice the radius

- **Statement.** In a Delaunay triangulation, every vertex within geodesic distance `r` of a vertex
  `m` is within distance `2r` of `m` along the edge graph (stretch factor below 2). Alg. 10 removes
  the inserted interior vertices within Dijkstra distance `ℓ_lm` of the midpoint `m` of a split
  boundary edge `lm`, which contains all those within geodesic distance `ℓ_lm/2` that Chew's
  algorithm must remove (IC paper Sec. 4.2; SPEC §5 item 4).
- **Reference.** G. Xia, *The Stretch Factor of the Delaunay Triangulation Is Less than 1.998*,
  SIAM J. Comput. 42(4), 1620–1659, 2013, <https://doi.org/10.1137/110832458>, Corollary 1, as
  cited by the IC paper Sec. 4.2, which states that Xia's planar proof (by triangle strips) "applies
  without modification" to intrinsic Delaunay triangulations of surfaces.
- **Source check.** Read: IC paper Sec. 4.2. Bibliographic data from Crossref. Not read: Xia's paper
  (paywalled); the corollary number and the transfer to surfaces are the IC paper's.
- **Why not formalized.** The stretch-factor bound is a long planar geometric argument; its
  extension to surfaces is asserted, not proved, in the IC paper.
- **What depends on it.** Only that refinement with boundary removes enough vertices; removing more
  is harmless (App. C.1). Not a correctness condition of the integer data.

### [formula: ic-mollify] Intrinsic mollification

- **Statement.** Adding the same `ε ≥ 0` to every edge length, `ε` the smallest value for which
  every triangle satisfies the triangle inequality with slack `δ` (IC: `δ = 10⁻⁵ h̄`), makes every
  face nondegenerate and changes the geometry only where a triangle was within `δ` of degenerate
  (IC paper Sec. 3.9; SPEC §7, I13).
- **Reference.** N. Sharp, K. Crane, *A Laplacian for Nonmanifold Triangle Meshes*, Computer
  Graphics Forum 39(5) (SGP 2020), 69–80, <https://doi.org/10.1111/cgf.14069>, intrinsic
  mollification; IC paper Sec. 3.9.
- **Source check.** Read: IC paper Sec. 3.9. Bibliographic data from Crossref. Not re-read: Sharp and
  Crane 2020b.
- **Why not formalized.** The formula itself **is** machine-checked, under its own id: the code's `ε`
  is the least feasible value and the mollified faces satisfy the strict triangle inequalities
  (`[formula: intrinsic-mollification]`, `DdgProofs/Bff/Quality.lean`; the Rust function is
  `ddg_core::mollify::mollify_edge_lengths`). What is cited is only the claim that this preprocessing
  makes the later floating-point operations robust, an empirical statement (IC paper Sec. 3.9,
  Table 1).
- **What depends on it.** Robustness on near-degenerate inputs, not correctness: the operation is
  off by default in the library and reported when on (SPEC §7).

### [formula: ic-vertex-removal-delta] Removing an inserted vertex of a Δ-complex by flips (open)

- **Statement (conjecture).** In a Δ-complex (self-edges, repeated vertices and triangles allowed),
  an inserted interior vertex of angle sum `2π` and degree `> 3` can be brought to degree 3 by flips
  of convex edges at it.
- **Reference.** IC paper App. D, the remark after Thm. D.1 ("we conjecture that an analogous
  technique could be applied"), pointing to N. Sharp, K. Crane, *You Can Find Geodesic Paths in
  Triangle Meshes by Just Flipping Edges*, ACM TOG 39(6), Article 249, 2020,
  <https://doi.org/10.1145/3414685.3417839>, App. A, Thm. A.5–A.8 (flips inside a wedge of angle
  `< π` in a Δ-complex terminate: either the wedge degree drops or the first twice-incident edge moves
  towards the start).
- **Source check.** Read: IC paper App. D; FlipOut App. A (Lemmas A.1–A.4, Thm. A.5–A.8). Crossref
  for the FlipOut bibliographic data.
- **Why not formalized: open, and the bounded attempt stopped here.** The angle half of Thm. D.1 is
  proved for any cyclic sequence of triangle corners at the vertex, which includes Δ-complex stars
  (`FlatStar.exists_flippable`, `[formula: ic-vertex-removal]`): some spoke has a convex
  quadrilateral. What fails in a Δ-complex is the combinatorial half: the spoke may be a loop (both
  ends at the vertex, twice incident), its neighbours `j_{k−1}`, `j_{k+1}` may be the vertex itself
  (the new edge is a loop at it), and `j_k` may have degree one, so the flip may be impossible or may
  not lower the degree. FlipOut's progress measure needs the wedge angle `< π`, which makes every
  inner corner `< π` for free; in a full `2π` star the inner corner of the chosen spoke is only
  known to be `≤ π` when it is chosen by the angle argument, and the first-twice-incident-edge
  potential is not shown to decrease under that choice. No counterexample was found either.
  WP9 (tier C) did not reopen it: the complex model now has the flip (`Cx.flip_valid`), so the
  combinatorial half is statable, but the angle half would have to be carried through a sequence of
  flips whose geometry the model does not have. The runs of SPEC §4.6 (after I27's merge) leave
  stuck removals only at vertices with an edge from the vertex to itself (70 of 4 609 in
  `removals_match_the_signposts_exactly`; 74 after Delaunay flipping), where neither a
  degree-lowering flip nor a merge applies; whether some other flip sequence removes them (whether
  they are counterexamples) is not known.
- **What depends on it.** Removal in refinement with boundary (SPEC §4.6): the library flips by the
  rule of I8/I24, bounds the loop by `10 deg` flips and refuses with `RemovalStuck` (state restored)
  when no flip helps; refusals are counted in the evaluation.

---

## Intrinsic Simplification

"IEM paper" is H.-T. D. Liu, M. Gillespie, B. Chislett, N. Sharp, A. Jacobson, K. Crane, *Surface
Simplification using Intrinsic Error Metrics*, ACM TOG 42(4), Article 118, 2023
(<https://doi.org/10.1145/3592403>; `../IntrinsicErrorMetric/IntrinsicErrorMetric.pdf`; equation,
figure and section numbers are the paper's). "CETM" is B. Springborn, P. Schröder, U. Pinkall,
*Conformal Equivalence of Triangle Meshes*, ACM TOG 27(3), Article 77, 2008
(`../IntrinsicErrorMetric/papers/springborn2008_cetm.pdf`). The contract is
`crates/ddg-intrinsic-simplification/SPEC.md` ("SPEC", decisions `S1`–`S24`, theorem table §13).

**What the machine-checked part covers, and what it does not.** Machine-checked
(`DdgProofs/IntrinsicSimplification/*.lean`, README section "Intrinsic Simplification"): one
flattening (the derivative behind Eq. 4 for stars with one, two or three corners per face at the
vertex, monotonicity, limits and roots of the angle sum with extended angles, conservation of
curvature), one removal (degree-3 and boundary layouts and barycentric coordinates, the strict
boundary analogue of the flip lemma), the per-removal metric (Eq. 5–9, the 1-Wasserstein reading of
the memoryless cost, centres of mass, the transport between frames of a flat region), Eq. 10, the
wrap of angular coordinates, Sec. 7.3.1's angles, the anisotropic factor, the Galerkin operators of
the hierarchy, the distortion measures and the QEM quadrics. Not covered: the greedy algorithm
(queue, re-evaluation, stopping), the data structure, Springborn's flip fallback (its success is
not guaranteed by the paper either), the Δ-complex flip lemma (`ic-vertex-removal-delta`, open),
the circumcircle property of CETM's projective interpolation (CETM App. D, Proposition; checked
in floating point by the oracle's `iem-selftest`), the quality of the result, and the three
results below.

### [formula: iem-wasserstein-memory] The memory-based metric approximates the transport cost of the fine mass

- **Statement.** After many removals, `Σ_j m_j ‖t_j‖` (Eq. 6 with the error vectors of Eq. 7/8)
  approximates the 1-Wasserstein distance between the initial fine mass distribution and the
  current coarse one. The paper says the cost "is only approximate since the 1-Wasserstein
  distance to the center of mass is not in general equal to the distance to the original fine
  distribution—but it is usually quite close" (Sec. 5.2), and on surfaces the transported vector
  "merely approximates the center of mass" (Sec. 5.3).
- **Reference.** IEM paper Sec. 5.1–5.3, Fig. 9; G. Peyré, M. Cuturi, *Computational Optimal
  Transport*, Foundations and Trends in Machine Learning 11(5–6), 2019, Ch. 2 (cited there for
  the 1-Wasserstein distance).
- **Source check.** Read: IEM paper Sec. 5 and Sec. 9 (the first of the three "basic
  approximations"). Not re-read: Peyré–Cuturi (standard definitions only).
- **Why not formalized.** It is an empirical approximation claim with no stated bound. What is
  machine-checked: for **one** removal the memoryless cost **is** the 1-Wasserstein distance
  between the old and the new distribution (`wasserstein_eq`, `[formula: iem-memoryless-cost]`),
  and in the plane Eq. 7 keeps `x_j + t_j` **exactly** at the centre of mass of all the fine mass
  merged into `j` (`eq7_centerOfMass`, `centerOfMass_combine`, `eq8_eq_eq7`,
  `[formula: iem-center-of-mass-update]`); the gap between "distance to the centre of mass" and
  "transport distance" is what the paper accepts.
- **What depends on it.** The meaning of the priority (which vertex goes first), not correctness:
  every removal is valid whatever its cost.

### [formula: iem-cetm-existence] Conformally equivalent metrics with prescribed angle sums: convex energy, uniqueness

- **Statement.** CETM's energy `E(u)` is convex, extended to all of `ℝ^V` by the extended angles
  (`C¹`, gradient `Θ̃ − Θ(u)`, Hessian the cotan Laplacian), so a solution with prescribed angle
  sums is the unique minimizer up to the constants of the problem when it exists; triangle
  inequality violations at the minimizer are handled by flips (Sec. 3.2), and the boundary
  condition `u = 0` gives the least distorted flat metric (App. E).
- **Reference.** CETM Sec. 2–3 (the energy, Eq. 8–10), Sec. 3.1 "Convex Optimization", Sec. 3.2
  "Violation of the Triangle Inequality", Sec. 4.2 (fixed `u`), App. E (the natural boundary
  condition, a smooth theorem).
- **Source check.** Read: CETM Sec. 2–4, App. D, E.
- **Why not formalized.** The multi-vertex energy and its convexity on `ℝ^V` (built from Milnor's
  Lobachevsky function) are not set up here. The **one-vertex case that the IEM paper uses** is
  machine-checked directly, without the energy: with extended angles `Θ_i(u)` is continuous and
  non-increasing, tends to `deg·π` and `0`, so a root exists for every target in `(0, deg·π)`,
  is unique where some face is a strict triangle, and any bracket of S9 contains one
  (`[formula: iem-flatten-monotone]`); Eq. 4 is Newton's step for it, its denominator the
  (corrected) derivative (`[formula: iem-flatten-newton]`). App. E's optimality of `u = 0` on the
  boundary of the star (IEM paper Sec. 4.1) is the reason for the choice, not a correctness
  condition.
- **What depends on it.** Nothing in the code beyond the one-vertex statements above.

### [formula: iem-delaunay-quality] Corner angles of at least 30° after Delaunay refinement of the coarse mesh

- **Statement.** Intrinsic Delaunay refinement of the coarse mesh yields corner angles of at least
  `30°` (hence at most `120°`) when the mesh is closed and every input angle sum is at least
  `60°` (IEM paper Sec. 6.3, "Post-refinement", citing Gillespie et al. 2021a). This is
  `[formula: ic-refinement-termination]` (IC paper Thm. C.1), applied to the coarse mesh as an
  input (SPEC §6.2); SPEC's default bound is `25°`, inside the guarantee.
- **Reference.** IEM paper Sec. 6.3; IC paper App. C.2, Thm. C.1 (see that entry above).
- **Source check.** Read: IEM paper Sec. 6.3; IC paper App. C as recorded under
  `ic-refinement-termination`.
- **Why not formalized.** As for `ic-refinement-termination`.
- **What depends on it.** Element quality of the optional post-refinement (SPEC §6.2), reported as
  best effort when the hypotheses fail.

## Repulsive Surfaces

"RS paper" is C. Yu, C. Brakensiek, H. Schumacher, K. Crane, *Repulsive Surfaces*, ACM TOG 40(6),
Article 268, 2021 (<https://doi.org/10.1145/3478513.3480521>); section, equation and figure
numbers are the paper's. The contract is `crates/ddg-repulsive-surfaces/SPEC.md` ("SPEC",
decisions `D1`–`D31`, formula table §12). "N&W" is J. Nocedal, S. J. Wright, *Numerical
Optimization*, 2nd ed., Springer 2006.

**What the machine-checked part covers, and what it does not.** There is no certifiable minimizer:
the library runs a descent method on a nonconvex energy. Machine-checked
(`DdgProofs/RepulsiveSurfaces/*.lean`, namespace `DDG.RS`, 35 `rs-` ids): the exponent
bookkeeping; the face quantities and dual areas; the derivatives of face area, unit normal and
barycenter along any differentiable motion; the kernel (Eq. 16), its invariances, homogeneity and
projector form, and its gradient for `α > 1` everywhere off the diagonal, **including `ν = 0`**;
the discrete energy (Eq. 10), its invariance under translations, orthogonal maps and face flips,
its homogeneity `c^{4+α−β}`, its exact gradient along any regular motion (nondegenerate faces,
distinct barycenters) and the identities `Σ g = 0`, `Σ f × g = 0`, `Σ f·g = (4+α−β)Ê`; the
Barnes–Hut energy with the leaf rule (exact up to the lumping of admissible clusters under an
exact cover; exact at `θ = 0`), the aggregate adjoint of D4, the box-distance lower bound; `D_f`
as the P1 gradient (with a counterexample to the printed unit-edge form); the operators of
Eq. 11–13 assembled by increments, symmetric PSD with the constants in their kernel; the
symmetrisation of `B⁰`; the App. A kernel-matrix form for `U` and `V` and the null space of the
compressed operator; the projector average; the preconditioner's order bookkeeping and its
symmetry, positivity and range; the Schur block inverse, direction, Newton correction (least
norm, row scaling, the sign) and the nested saddle; the constraint Jacobians (area, volume with
the factor 3 on closed fans, per-component barycenter and its exact translation projection); the
penalties of §5.2 and the structure of the Willmore derivative; the reference's Armijo test; the
circumcenter weights, smoothing weights and Delaunay test of the remesher; conservative
advancement with moving triangles; the L-BFGS two-loop recursion and positive definiteness of the
update. **Not** machine-checked: the entries below; anything about the discrete energy as an
approximation of the smooth one; the derivatives of the cotan weights and masses inside the
Willmore gradient (finite differences, SPEC §13.2); the turning-angle gradient at `θ_v = 0`
(the formula is proved for `0 < θ_v < π`); the error of Barnes–Hut and of the block cluster tree;
the convergence of GMRES, of the Newton projection, of the flow and of the comparison
optimizers; the data structures (BVH, BCT, remeshing mesh) and every tolerance; rounding.

**SPEC errors found while proving (both fixed in the SPEC).** The first SPEC version printed the
Newton correction as "`h = A′⁻¹ C_Nᵀ S⁻¹ Φ(f)` subtracted (`C_N h = Φ`)". With
`S = −C_N A′⁻¹ C_Nᵀ` that `h` satisfies `C_N h = −Φ` (`RS.schur_correction`), so it must be
**added** (paper Eq. 21, "which we add to the updated embedding"); subtracting it doubles the
linearized violation (`RS.schur_correction_sign`). The SPEC (§6.2, since 5ead321) now defines the
Rust step as the negation, `h = −A′⁻¹ C_Nᵀ S⁻¹ Φ`, with `C_N h = Φ`, applied as `f ← f − h`: the
same update. SPEC §6.3 tagged `rs-fractional-preconditioner` without a §12 row; the row is added,
and the id is proved (`RS.aug_precond_props`, the default augmented form of D33; the printed form
`RS.frac_precond_props` is kept untagged as the reference mode).

### [formula: rs-exponent-range] The energy space of the tangent-point energy of surfaces

- **Statement.** For an `n`-dimensional set `Σ ⊂ R^m` with finite tangent-point energy
  `ℰ_p(Σ)`, `p > 2n`, `Σ` is an embedded `C^{1,1−2n/p}` submanifold (Strzelecki–von der Mosel);
  `ℰ_p(Σ) < ∞` iff `Σ` is an embedded submanifold of class `W^{s,p}`, `s = 2 − n/p` (Blatt).
  For surfaces (`n = 2`), `p > 4` makes the energy a barrier against self-contact, and the
  energy space is `W^{s,p}` with `s ∈ (3/2, 2)`.
- **Reference.** RS paper Sec. 2.1 (units `m^{2n−p}`, "`p > 2n` is required for the energy to be
  truly repulsive") and Sec. 2.3; P. Strzelecki, H. von der Mosel, *Tangent-point repulsive
  potentials for a class of non-smooth m-dimensional sets in R^n. Part I*, J. Geom. Anal. 23(3),
  1085–1139, 2013; S. Blatt, *The energy spaces of the tangent point energies*, J. Topol. Anal.
  5(3), 261–270, 2013. The two-exponent form `(α, β)`: S. Blatt, P. Reiter, *Regularity theory for
  tangent-point energies: the non-degenerate sub-critical case*, Adv. Calc. Var. 8(2), 93–116,
  2015 (cited in Sec. 2.1).
- **Source check.** Read in the RS paper (Sec. 2.1, 2.3); the cited articles were not re-read, and
  the paper gives no theorem numbers for them.
- **Why not formalized.** Fractional Sobolev–Slobodeckij spaces on surfaces and the regularity
  theory of the energy are not in Mathlib.
- **What depends on it.** Only the admissible exponent range and the "subcritical" flag of
  `Exponents` (SPEC §2.2, D3): `s = (β − 2)/α ∈ (1, 2)` (the bookkeeping is machine-checked,
  `rs-exponents`), subcritical `β − α ≤ 4`. The discrete formulas and their proofs hold for any
  real exponents.

### [formula: rs-differential-order] The differential has order `2s`; `B` matches it

- **Statement.** `dℰ_p` maps `W^{s,p}` to its dual `W^{−s,p′}`; in this sense it is a nonlocal
  differential operator of order `2s = 2(2 − n/p)` (`∈ (3, 4)` for surfaces and `p > 4`); the
  operator `B` of Eq. 7 with `σ = s − 1` has the same order `2s`, and `B⁰` (Eq. 8) lower order.
- **Reference.** RS paper Sec. 2.3 ("Order of the Differential") and Sec. 2.4, Eq. 5–8.
- **Source check.** Read in the RS paper.
- **Why not formalized.** A heuristic statement about smooth operators between fractional spaces;
  it motivates the choice of the metric and is not used by any formula.
- **What depends on it.** The efficiency of the `Hs` descent, not its correctness: the metric is
  symmetric PSD (`rs-sobolev-forms`), and every accepted step decreases the energy (checked by
  the line search, SPEC §7.2).

### [formula: rs-fractional-operator] The fractional Laplacian and `L_σ` have order `2σ`

- **Statement.** For `0 < 2σ < 2`, `⟨(−Δ)^σ u, v⟩ = c ∬_{R^n×R^n} (u(x) − u(y))(v(x) − v(y))
  |x − y|^{−(2σ+n)} dx dy` (up to a constant), and its analogue `L_σ` (Eq. 6) on the immersed
  surface is an operator of order `2σ`.
- **Reference.** RS paper Sec. 2.4, Eq. 5 and 6, citing M. Kwaśnicki, *Ten equivalent definitions
  of the fractional Laplace operator*, Fract. Calc. Appl. Anal. 20(1), 7–51, 2017.
- **Source check.** Read in the RS paper; Kwaśnicki's article was not re-read (journal data from
  memory, not from the paper's reference list, which is cut off in the text extraction).
- **Why not formalized.** Fractional Laplacians and their integral representations are not in
  Mathlib.
- **What depends on it.** The interpretation of `L_{2−s}` in the preconditioner and of the order of
  `B`. The discrete operators are defined by their formulas (Eq. 11–13), whose quadratic forms,
  symmetry, positivity and kernels are machine-checked (`rs-sobolev-forms`,
  `rs-kernel-matrix-form`).

### [formula: rs-preconditioner-approx] `(−Δ)⁻¹ L_{2−s} (−Δ)⁻¹` approximates the inverse of the high-order term

- **Statement.** Since `B` has order `2s` and `L_{2−s}` order `2(2 − s)`, the operator
  `Ã⁻¹ = (−Δ_M)⁻¹ L_{2−s} (−Δ_M)⁻¹` has order `−2s` and approximates `B⁻¹` well enough to be an
  effective preconditioner for GMRES on `A` (plus auxiliary terms), though not by itself a good
  inner product.
- **Reference.** RS paper Sec. 5.2 (the factorization `(−Δ)^{−s} = (−Δ)⁻¹(−Δ)^{2−s}(−Δ)⁻¹` and
  the replacement of `(−Δ)^{2−s}` by `L_{2−s}`, "as per Sec. 2.4"), and Fig. 13 for the measured
  iteration counts.
- **Source check.** Read in the RS paper.
- **Why not formalized.** A statement about spectral equivalence of nonlocal operators under
  refinement; the paper gives a heuristic and measurements. The exact matrix identity
  `W⁻¹ (M Δ^{2−s}) W⁻¹ = Δ^{−s} M⁻¹` is machine-checked (`rs-preconditioner-order`), as are the
  symmetry and nonsingularity of the augmented preconditioner and its reproduction of the constraint
  rows (`rs-fractional-preconditioner`, `RS.aug_precond_props`; positivity of the printed form is
  `RS.frac_precond_props`).
- **What depends on it.** The number of GMRES iterations (SPEC §6.3, §14.5), not the solution:
  GMRES stops on the true residual of `A′` (right preconditioning, D13).

### [formula: rs-preconditioner-scale] The constant between `W⁻¹ L_{2−s} W⁻¹` and `B⁻¹`

- **Statement.** On a plane, for `u = e^{ik·x}`, the quadratic forms of `B` (Eq. 12, kernel
  `|x − y|^{−2s}`, `σ = s − 1`), of the vertex-pair `L_{2−s}` (kernel `|x − y|^{−2(2−s)−2}`) and of the
  cotan stiffness `W` have the symbols `2 I(σ) |k|^{2s}`, `2 I(2−s) |k|^{2(2−s)}` and `|k|²` per unit
  area, with `I(t) = ∫_{R²} (1 − cos z₁) |z|^{−2−2t} dz = π Γ(1−t) / (t 4^t Γ(1+t))`, `t ∈ (0, 1)`;
  hence `W⁻¹ L_{2−s} W⁻¹ B` has the symbol `4 I(σ) I(2−s)` (≈ 199.86 at `s = 5/3`), and the augmented
  preconditioner divides `L_{2−s}` by `w_TP · 4 I(σ) I(2−s)` (SPEC §6.3, `FractionalScale::Symbol`).
- **Reference.** E. Di Nezza, G. Palatucci, E. Valdinoci, *Hitchhiker's guide to the fractional
  Sobolev spaces*, Bull. Sci. Math. 136 (2012), Sec. 3 (Prop. 3.3: the Gagliardo seminorm equals
  `2 C(n,s)⁻¹ ∫ |ξ|^{2s} |û|²` with `C(n,s)⁻¹ = ∫ (1 − cos ζ₁) |ζ|^{−n−2s} dζ`, and the closed form
  of `C(n,s)`), with `n = 2`.
- **Source check.** Not re-read in this session; the closed form of `I(t)` was checked against an
  adaptive quadrature at `t = 1/3, 2/3, 0.4286, 0.5714` (agreement to 6 digits;
  `crates/ddg-repulsive-surfaces/NOTES-SOLVERFIX.md`).
- **Why not formalized.** A Fourier-symbol computation for the continuous operators; the discrete
  operators only approximate it (measured: the preconditioned spectrum of a 900-face torus lies in
  0.59–1.94 with it, in 120–300 without).
- **What depends on it.** Only the iteration count (the relative position of the constraint part,
  eigenvalue 1, and the fractional part of the preconditioned spectrum), never the solution.

### [formula: rs-aqp] Accelerated quadratic proxy

- **Statement.** AQP minimizes a quadratic proxy (here the cotan Laplacian `W`) with a Nesterov
  acceleration step from the previous two configurations and a line search on the true energy;
  it supports linear constraints only.
- **Reference.** S. Z. Kovalsky, M. Galun, Y. Lipman, *Accelerated Quadratic Proxy for Geometric
  Optimization*, ACM TOG 35(4), Article 134, 2016, Sec. 2 (as cited in RS paper Sec. 8.2.1).
- **Source check.** RS paper Sec. 8.2.1 read; Kovalsky et al. not re-read (the momentum constants
  are the SPEC's D28 reading of it).
- **Why not formalized.** A convergence-acceleration heuristic for nonconvex energies; there is no
  formula whose correctness the comparison depends on beyond the saddle solve (machine-checked,
  `rs-schur-direction`, `rs-nested-saddle`).
- **What depends on it.** Only the comparison runs of SPEC §14.2 (D28); not the main flow.

### [formula: rs-bcqn] Blended cured quasi-Newton

- **Statement.** BCQN interpolates between the `H¹` (Laplacian) step and the `H¹`-initialized
  L-BFGS step with a weight driven by the measured curvature, and adds inversion barriers
  (omitted here, RS paper Sec. 8.2.1).
- **Reference.** Y. Zhu, R. Bridson, D. M. Kaufman, *Blended Cured Quasi-Newton for Distortion
  Optimization*, ACM TOG 37(4), Article 40, 2018 (as cited in RS paper Sec. 8.2.1).
- **Source check.** RS paper Sec. 8.2.1 read; Zhu et al. not re-read (the blend with a
  spectral-norm estimate is the SPEC's D27 reading of it).
- **Why not formalized.** An optimization heuristic; the L-BFGS part is machine-checked
  (`rs-lbfgs-two-loop`).
- **What depends on it.** Only the comparison runs of SPEC §14.2 (D27).

### [formula: rs-energy-consistency] The discrete energy approximates the smooth one

- **Statement.** For a smooth embedded surface `f` and sampled triangulations `f_h`, the discrete
  energy `Ê(f_h)` of Eq. 10 converges to `ℰ_p(f)` at an empirical rate between `O(h)` and
  `O(h²)`, and the Barnes–Hut energy `Ẽ` converges to `Ê` as `O(θ²)`.
- **Reference.** RS paper Sec. 8.1 and Fig. 10 (and the supplemental examples cited there).
- **Source check.** Read in the RS paper.
- **Why not formalized.** An empirical observation; the paper gives no proof (it notes that face
  normals alone would suggest `O(h)`).
- **What depends on it.** That minimizers of the discrete energy resemble those of the smooth one;
  SPEC §14.1 re-measures it. The flow minimizes the discrete energy as defined (`rs-energy`), whose
  gradient is machine-checked (`rs-energy-gradient`); at `θ = 0` the Barnes–Hut energy is exact
  (`rs-bh-energy`).

### [formula: rs-chen-holst] The remeshing scheme

- **Statement.** Splitting long and collapsing short edges, Delaunay flips and tangential
  smoothing toward the area-weighted circumcenters (optimal Delaunay triangulation smoothing)
  improve triangle quality while keeping the surface.
- **Reference.** RS paper Sec. 6, citing L. Chen, M. Holst, *Efficient mesh optimization schemes
  based on Optimal Delaunay Triangulations*, Comput. Methods Appl. Mech. Engrg. 200(9), 967–984,
  2011.
- **Source check.** RS paper Sec. 6 read; Chen–Holst not re-read.
- **Why not formalized.** Mesh quality improvement is a heuristic; its local formulas are
  machine-checked (`rs-circumcenter`, `rs-smoothing-weights`, `rs-delaunay-test`).
- **What depends on it.** The triangle quality during long flows (SPEC §10); every remeshing result
  is revalidated as a manifold `Surface`, and foldovers are refused by explicit checks.

---

## Points as Tori

"PAT paper" is N. Feng, I. Gkioulekas, K. Crane, *Points as Tori: Fast Pointwise Signed Distance
for Point Clouds*, ACM TOG 45(4), Article 53, 2026 (<https://doi.org/10.1145/3811385>); section,
equation and figure numbers are the paper's. The contract is `crates/ddg-points-as-tori/SPEC.md`
("SPEC", formula table §13); "N§k" is the reference survey `crates/ddg-points-as-tori/NOTES-REFERENCE.md`.

**What the machine-checked part covers, and what it does not.** Machine-checked
(`DdgProofs/PointsAsTori/*.lean`, namespace `DDG.PAT`): the neighborhood scale (the median
and its homogeneity), the reference frame rule on both branches, the translation and scale
invariance of the features; the Monge patch, its fundamental forms, `H` and `K` of Eq. 19, the
eigenvector pair of the reference and Eq. 20, the Rayleigh quotients, the principal directions and
Eq. 21; the radii, the sign rule Eq. 22, the axis, and the scale covariance of the whole fit; the
torus SDF as the Euclidean distance to the torus (for `R ≥ |r|`; a lower bound for spindle tori),
its gradient and Laplacian off the axis and the core circle, the stable form; the contact of the
fitted torus with the patch (value, normal and the torus's principal curvatures, including the
cases where the paper's construction fails, below); the blend's shift invariance, weight range,
convexity, frozen and exact gradients (App. D) and its Laplacian; and the further ids listed in
SPEC §13 as proved, in the files named there. **Not** machine-checked: the entries below; the
network and everything learned; neighbor search; the clamps of D9 beyond the statements in
`Fitting.lean`; the one-sided rules on the singular sets (D12); floating-point rounding.

**Findings while proving (reported to the SPEC owner, SPEC not edited here).**
1. *Saddle tori are mirrored.* With the paper's center `c = Q*(0,0) − sign(T)/|κmin| n*(0,0)` and
   Eq. 22, the torus's principal curvatures at `o = Q*(0,0)` w.r.t. `n*` are `(κmax, κmin)` only
   for `κmin κmax > 0`; for saddles they are `(−κmax, −κmin)`, and the SDF gradient at `o` is
   `−n*` (`PAT.torus_contact_curvatures`, `PAT.torus_contact_normal`). The center
   `c = o + n*/κmin` (with the same radii and Eq. 22 for the SDF sign) osculates in all cases with
   `∇g(o) = +n*`; both rules together are `PAT.fit_contact` (`pat-fit-contact`: `Geometric` gives
   `∂g/∂n*(o) = +1`, `Reference` gives `σ = −1` exactly at saddles). The reference code
   (`torus_distance.cpp:781-807`) is the paper's, and the pretrained network was trained with it
   (SPEC §18: `SignRule::Geometric` is the default since 2026-11-04; `SignRule::Reference` reproduces
   the authors' program and training).
2. *Degenerate eigenvector pair.* `(v1 ± √Δ, v2)` vanishes for one sign whenever `v2 = 0`,
   `v1 ≠ 0`, e.g. `a10 = a01 = a11 = 0`, `a20 ≠ a02` (`PAT.principal_pair_degenerate`); the
   reference's `1e-10` clamp then returns the other eigenvector twice.
3. *Contact needs `κ₊κ₋ ≠ 0`.* With `sgn(0) = 0` (D9), an exact zero curvature puts `o` on the
   core circle: `g(o) = −sgn(τ r)|r| ≠ 0` (`PAT.torus_contact_fails_sign_zero`).
4. *Weight range.* The weights are in `[e^{−C}, e^{C}]` and `W ≥ |𝒩| e^{−C}`
   (`PAT.blend_weight_range`, `PAT.blend_weight_sum_pos`); the nearest point's weight is `≥ 1` only
   when `ρ_min ≤ σ`, not in general as SPEC §7.4 says.

### [formula: pat-network] The network architecture

- **Statement.** The coefficient network is the stack of SPEC §6.1: per-point MLP embedding, eight
  pre-norm transformer blocks (multi-head self-attention and an MLP, residual connections, layer
  normalization) and a linear read-out of the first token (PAT paper Sec. 4.3, Fig. 15; N§2).
- **Reference.** A. Vaswani et al., *Attention is all you need*, NeurIPS 2017; J. L. Ba,
  J. R. Kiros, G. E. Hinton, *Layer normalization*, arXiv:1607.06450, 2016; the Flax
  `nnx.MultiHeadAttention` and `nnx.LayerNorm` modules used by the reference (N§2.3).
- **Source check.** PAT paper Sec. 4.3 and N§2 read; the cited articles were not re-read.
- **Why not formalized.** It is a definition of a learned function, not a mathematical claim; its
  output has no property beyond being six real numbers. Correctness of the port is a parity
  question (SPEC §15.1 item 4).
- **What depends on it.** The coefficients `a`. Every proved statement downstream holds for any
  coefficients.

### [formula: pat-laplace-method] Laplace's method and the self-normalized limit

- **Statement.** For continuous `h > 0` and twice-differentiable `φ` with a unique minimizer `x*` on
  `Ω`, `∫_Ω h e^{−λφ} ~ (2π/λ)^{d/2} det(∇²φ(x*))^{−1/2} h(x*) e^{−λφ(x*)}` as `λ → ∞` (Eq. 15), so
  `−(1/λ) log ∫ h e^{−λ|x − z|} dA(z) → min_z |x − z|` (Eq. 16) and the self-normalized average
  `∫ g e^{−λ|x−z|} / ∫ e^{−λ|x−z|}` (Eq. 17) tends to `g` at the closest point.
- **Reference.** PAT paper Sec. 2.2, Eq. 15–17; L. C. Evans, *Partial Differential Equations*, AMS
  1998, §4.5; C. M. Bender, S. A. Orszag, *Advanced Mathematical Methods for Scientists and
  Engineers I*, Springer 1999, §6.4.
- **Source check.** PAT paper Sec. 2.2 read; Evans and Bender–Orszag not re-read.
- **Why not formalized.** Asymptotic analysis of integrals over surfaces; Mathlib has no Laplace
  method for surface integrals. The library evaluates the discrete sum (Eq. 25) and needs only its
  exact algebra (`pat-blend-*`, proved).
- **What depends on it.** Only the motivation for Eq. 25 and the choice of large `λ`.

### [formula: pat-signed-hopf-cole] The signed Hopf–Cole transformation

- **Statement.** The screened Poisson problems of Sec. 2 (Eq. 8–14) and their boundary-integral
  forms give, through `u = −(1/λ) log w` (signed), an approximation of the signed distance that
  converges as `λ → ∞`; App. A derives the boundary conditions and the link to winding numbers,
  App. B the convergence.
- **Reference.** PAT paper Sec. 2.1, Eq. 8–14, App. A–B; Y. Lipman, *Phase transitions, distance
  functions, and implicit neural representations*, ICML 2021 (arXiv:2106.07689); S. R. S.
  Varadhan, *On the behavior of the fundamental solution of the heat equation with variable
  coefficients*, Comm. Pure Appl. Math. 20, 1967.
- **Source check.** PAT paper Sec. 2.1 and the statements of App. A–B read; Lipman and Varadhan not
  re-read.
- **Why not formalized.** PDE theory on domains with boundary (potentials, jump conditions,
  viscosity limits) that Mathlib does not have.
- **What depends on it.** Nothing computed: the theory motivates the method. SHC (Eq. 28) is the
  discretization of Eq. 14 and its algebra is proved (`pat-shc`).

### [formula: pat-winding-dipole] The point-cloud winding number as a sum of dipoles

- **Statement.** For an oriented point cloud with areas `A_i`, the generalized winding number is
  `w(x) = Σ_i A_i ⟨n_i, p_i − x⟩ / (4π |p_i − x|³)`, the discretization of the solid-angle
  (double-layer) integral; `w ≈ 1` inside and `≈ 0` outside a closed well-sampled surface.
- **Reference.** G. Barill, N. Dickson, R. Schmidt, D. I. W. Levin, A. Jacobson, *Fast winding
  numbers for soups and clouds*, ACM TOG 37(4), 2018, Sec. 3, Eq. 6; A. Jacobson, L. Kavan,
  O. Sorkine-Hornung, *Robust inside-outside segmentation using generalized winding numbers*,
  ACM TOG 32(4), 2013.
- **Source check.** Not re-read (the formula as in SPEC §9 and N§7.6).
- **Why not formalized.** The approximation property (a quadrature of the solid angle) is
  analytic; the per-term kernel is the dipole `−n·∇G` of `G = 1/(4π|y|)`, whose derivatives are
  proved under `pat-fwn-expansion`.
- **What depends on it.** The GWN baseline and the sign of the mesh ground truth (SPEC §9).

### [formula: pat-point-area] Point areas of Barill et al.

- **Statement.** The area of `p_i` is the area of its Voronoi cell in the Delaunay triangulation of
  its `k` nearest neighbors projected to the tangent plane.
- **Reference.** Barill et al. 2018 (above), Sec. 5.
- **Source check.** Not re-read.
- **Why not formalized.** A heuristic quadrature weight; no exact property is claimed.
- **What depends on it.** The area weights of GWN, FWN, SHC and SSPD (baselines only).

### [formula: pat-solid-angle] The solid angle of a triangle

- **Statement.** For corners `a, b, c` relative to the viewpoint, the solid angle is
  `Ω = 2 atan2(det[a, b, c], |a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b|)`; the generalized winding
  number of a mesh at `x` is `Σ Ω / 4π`.
- **Reference.** A. van Oosterom, J. Strackee, *The solid angle of a plane triangle*, IEEE Trans.
  Biomed. Eng. BME-30(2), 1983, 125–126, Eq. (8); the winding number: A. Jacobson, L. Kavan,
  O. Sorkine-Hornung, *Robust inside-outside segmentation using generalized winding numbers*,
  ACM TOG 32(4), 2013, Eq. 3.
- **Source check.** Not re-read.
- **Why not formalized.** Mathlib has no solid angles of triangles (cf. `spherical-triangle-area`).
- **What depends on it.** The sign of the mesh ground truth (exact winding number) and the
  triangle elements of the fast winding number (SPEC §9).

### [formula: pat-triangle-sampling] Uniform sampling of a triangle

- **Statement.** For independent `u, v ~ U(0,1)`, the point `√u(1−v) a + √u v b + (1 − √u) c` is
  uniformly distributed on the triangle `abc`; a face is drawn with probability proportional to
  its area.
- **Reference.** R. Osada, T. Funkhouser, B. Chazelle, D. Dobkin, *Shape distributions*, ACM TOG
  21(4), 2002, Sec. 4.2 (`P = (1−√r₁)A + √r₁(1−r₂)B + √r₁ r₂ C`; the corners are relabeled here,
  as in the authors' sampler).
- **Source check.** Not re-read; the barycentrics are those of the authors' mesh sampler
  (`NOTES-REFERENCE.md` §9, `MS`/`S3.cpp`).
- **Why not formalized.** A statement about a random variable (a change of variables on the unit square);
  nothing downstream depends on it beyond the training distribution.
- **What depends on it.** The training queries (`train::queries`); no inference result.

### [formula: pat-query-sampling] Training queries of a neighborhood

- **Statement.** For one training row, the queries are drawn from the faces whose bounding boxes
  meet the box of the row's neighborhood (the box of the neighbor positions, widened to a minimum
  extent), by area, and perturbed as the options say (a definition of the training distribution
  as coded, SPEC §10.3).
- **Reference.** The authors' training-data sampler (`points-as-tori`, MIT)
  as recorded in `NOTES-REFERENCE.md` §9 (`MS` = `training/mesh_sampling_jax.py`) and `NOTES-WP7.md` item 2; the paper's Sec. 4.3 and
  App. C describe the distribution in words.
- **Source check.** Reference code read (NOTES-WP7.md); paper text not re-read.
- **Why not formalized.** A definition of a sampling procedure, not a mathematical property; it
  has no formula to prove.
- **What depends on it.** The training distribution only.

### [formula: pat-regularized-wn] The regularized winding number

- **Statement.** The dipole kernel multiplied by `S(t) = erf t − (2t/√π) e^{−t²}`, `t = ρ/ε`, is
  smooth at `ρ = 0` and equals the dipole kernel up to `O(e^{−t²})` (the field of a Gaussian
  smoothed dipole).
- **Reference.** H. Chen, B. Miller, I. Gkioulekas, *3D reconstruction with fast dipole sums*,
  ACM TOG 43(6), Article 192, 2024 (PAT paper "Chen et al. 2024a").
- **Source check.** Not re-read.
- **Why not formalized.** A smoothing heuristic for a baseline; the statement is analytic.
- **What depends on it.** The regularized-WN baseline only.

### [formula: pat-lse] LogSumExp distance

- **Statement.** `d_λ(x) = −(1/λ) log Σ_i e^{−λ|x − p_i|}` satisfies
  `min_i |x − p_i| − (log n)/λ ≤ d_λ(x) ≤ min_i |x − p_i|`, a smooth unsigned distance.
- **Reference.** A. Madan, D. I. W. Levin, *Fast evaluation of smooth distance constraints on
  co-dimensional geometry*, ACM TOG 41(4), 2022; P. Blanchard, D. J. Higham, N. J. Higham,
  *Accurately computing the log-sum-exp and softmax functions*, IMA J. Numer. Anal. 41(4), 2021
  (the shift).
- **Source check.** Not re-read.
- **Why not formalized.** Baseline only; the bound is elementary but nothing in the library relies
  on it. The shift invariance is the same algebra as `pat-blend-shift` (proved).
- **What depends on it.** The LSE baseline only.

### [formula: pat-adam] Adam

- **Statement.** `m ← β₁m + (1−β₁)g`, `v ← β₂v + (1−β₂)g²`, `θ ← θ − η m̂/(√v̂ + ε)` with bias
  corrections `m̂ = m/(1−β₁ᵗ)`, `v̂ = v/(1−β₂ᵗ)`.
- **Reference.** D. P. Kingma, J. Ba, *Adam: a method for stochastic optimization*, ICLR 2015
  (arXiv:1412.6980).
- **Source check.** Not re-read; the update is the framework's.
- **Why not formalized.** An optimizer heuristic with no convergence claim used here.
- **What depends on it.** Training (SPEC §10) only.

### [formula: pat-offset] Offsets and constructive solid geometry of SDFs

- **Statement.** For an exact SDF `φ` of a set `Ω`, the level set `{φ = c}` is the boundary of the
  `c`-offset of `Ω` (`c > 0`) or of its erosion (`c < 0`), and `φ − c` is its exact SDF on the side
  away from the medial axis. `min(φ_A, φ_B)` is the exact SDF of `A ∪ B` outside the union and a
  bound inside; `max` likewise for intersections (Ricci). The approximate `φ` of Eq. 25 inherits
  these only approximately.
- **Reference.** A. Ricci, *A constructive geometry for computer graphics*, The Computer Journal
  16(2), 1973; I. Quilez, *distance functions*, <https://iquilezles.org/articles/distfunctions/>.
- **Source check.** Not re-read.
- **Why not formalized.** The exact-SDF statements concern sets and their metric projections; the
  library applies them to an approximate field. The elementary part is proved under
  `pat-boolean` (the sign rules and the preservation of 1-Lipschitz bounds by `min`/`max`).
- **What depends on it.** The offsets, Booleans and morphology applications (SPEC §8.1–8.3), which
  report the field they compute and claim no exactness.

### [formula: pat-quilez-2d] The 2D analytic SDFs

- **Statement.** The sixteen 2D SDFs of SPEC §11 (circle, pie, arc, segment, vesica, box, cross,
  isosceles triangle, triangle, hexagon, ellipse, heart, pentagon, moon, trapezoid, quad) are the
  exact signed distances to their shapes, with the stated gradients.
- **Reference.** I. Quilez, *2D distance and gradient functions*, 2019,
  <https://iquilezles.org/articles/distgradfunctions2d/>; PAT paper App. E.
- **Source check.** Not re-read.
- **Why not formalized.** Sixteen case analyses of planar nearest points; they are test data, and
  the Rust tests check them against brute-force distances to dense boundary samples.
- **What depends on it.** The analytic test shapes (SPEC §11) and their ground truth. The lifting
  to 3D (revolution, extrusion) and its gradients are proved (`pat-revolution`,
  `pat-extrusion`, `pat-extrusion-gradient`); their exactness given exact 2D SDFs is cited here
  (Quilez, *distance functions*, above).

### Not formalized for proved ids (cross-references)

- `pat-fwn-expansion`: the derivatives of `G(y) = 1/(4π|y|)` used by the order-0/1/2 terms are
  proved; the error bound of the truncated expansion for `β`-admissible clusters is Barill et al.
  2018, Sec. 4 (not re-read), and the octree is not covered.
- `pat-torus-sdf`: for `R < |r|` (spindle tori) only the lower bound `|g(x)| ≤ dist(x, T)` is
  proved; Eq. 23 is then the distance to the tube of the generating circle (Quilez 2025,
  *3D distance and gradient functions*).

---

## Out of scope of this package

The end-to-end pipeline (`crates/ddg-uv`: charts, packing), the Python bindings and the Blender
add-on are not covered. For `crates/ddg-bff` and `crates/ddg-optimal-cones` the formulas listed in
the README are covered; their control flow (iteration limits, safeguards, fallbacks, tolerances)
is not. For `crates/ddg-variational-cuts` the same holds: the formulas of one flow step are
covered, the extraction of the cut and of the cut mesh (which faces are split, the choice of the
quad diagonal, the bookkeeping of regions and components), the run loop around the stop rule,
the `N`-region variant (SPEC §8) as a whole and the helpers for initial level sets are not. For
`crates/ddg-intrinsic` the closed forms of the sixteen proved `it-` ids are covered; the data structure,
the algorithms built on them (flipping, tracing, insertion, overlay, refinement, export) and
their tolerances are not. For `crates/ddg-developability`
the formulas of the thirty proved `dev-` ids are covered; the flow driver, the remeshing and cutting
algorithms, the evaluation loop and every tolerance are not. For `crates/ddg-repulsive-curves` the
formulas of the forty-four proved `rc-` ids (and the six `rc-near-` ids of WP16) are covered; the session and step driver, the line
searches as algorithms, the BVH, block cluster tree and multigrid data structures, the scene
reader, the evaluation and every tolerance are not. For `crates/ddg-integer-coords` the integer
formulas, one face (the face model) and the local statements of the section "Integer Coordinates"
are covered, and the tier C statements of the complex model (trace termination, Alg. 3's index,
the flip with Eq. 3 and 4, the merge); the data structure, the walks of Alg. 2–4 as programs, the
common subdivision, refinement and every floating-point decision are not.

### [formula: pat-render-ggx] The GGX specular term of the demo shader

- **Statement.** `D F V` with `D = α²/(π((n·h)²(α² − 1) + 1)²)` (GGX / Trowbridge–Reitz
  normal distribution), `F = F0 + (1 − F0)(1 − l·h)⁵` (Schlick), and the Smith visibility with
  Schlick's approximation `V = Π 1/((n·w)(1 − k) + k)`, `k = α/2`, `α = roughness²`.
- **Reference.** B. Walter et al., *Microfacet Models for Refraction through Rough Surfaces*, EGSR
  2007; C. Schlick, *An Inexpensive BRDF Model for Physically-based Rendering*, 1994; the shader
  of the authors' demo (`demo/shaders_3d.py` of `points-as-tori`, MIT, from shadertoy `MlB3DV`).
- **Source check.** The constants (metallic 0.07, roughness 0.1, `F0 = 0.2`, factors 0.3 and 0.6)
  were read from the demo's `normalShade` and `shade`.
- **Why not formalized.** A shading model chosen for appearance; it has no claim to prove.
- **What depends on it.** `render::shade` only (the colour of the images of Fig. 12).

### [formula: pat-render-fresnel] The rim and Fresnel weights of the demo shader

- **Statement.** `(1 − |ω·n|)^80` (the rim term of `fresnelShade`) and `(1 − n·v)₊^p`,
  `p = mix(5, 3.5, metallic)` (the reflection weight of `normalShade`).
- **Reference.** The authors' demo shader (as above); the exponent form is Schlick's Fresnel.
- **Source check.** Read from `fresnelShade` and `normalShade`.
- **Why not formalized.** As `pat-render-ggx`.
- **What depends on it.** `render::shade` only.
