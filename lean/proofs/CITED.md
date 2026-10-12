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

### [formula: cotan-fem-stiffness] The cotan matrix as a finite element stiffness matrix

- **Statement.** `L` is the stiffness matrix of the Laplace–Beltrami operator for piecewise
  linear elements, `L u = M f` discretizes `Δu = f`, and solutions converge under refinement.
- **Reference.** Paper Sec. 5.1.3; MacNeal 1949; Pinkall, Polthier 1993.
- **Why not formalized.** Proved here: on each triangle `Area·|∇f|² = ½ Σ cot θ_k (f_i − f_j)²`
  for affine `f` (`[formula: cotan-formula]`) and `uᵀLu = Σ_faces` of these
  (`[formula: cotan-laplacian]`). Not formalized: the identification of `Area·|∇f|²` with the
  Lebesgue integral `∫_T |∇f|² dA` (the integrand is constant), the lumped mass matrix as a
  quadrature of the Galerkin mass matrix, and all convergence statements.
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
right-hand side, the Euler characteristic after cutting, and the definitions of the distortion
measures. Cited:

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

### [formula: bff-conformal-energy] Conformal energy and harmonic conjugation

- **Statement.** The least-squares conformal energy of `f = a + ib` is the Dirichlet energy
  minus the signed area of the image, `E_C(a, b) = ½(aᵀAa + bᵀAb) − 𝒜(a, b)` (BFF Eq. 15 is twice
  this, with `aᵀUb = −𝒜(a, b)` of its Eq. 16); it vanishes exactly
  for (discrete) conformal maps, and for fixed `a` its minimizer `b` is the "most conjugate"
  harmonic function.
- **Reference.** BFF Sec. 4.3.3, Eqs. 15–17; Hutchinson, *Computing conformal maps and minimal
  surfaces* (1991); Lévy, Petitjean, Ray, Maillot, *Least squares conformal maps*, ACM TOG 21(3),
  2002; Mullen, Tong, Alliez, Desbrun, *Spectral conformal parameterization*, SGP 2008, Sec. 2.2.
- **Why not formalized.** The identification of the discrete energy with `∫|∂̄f|²` and the
  statement about conformal maps are smooth/finite-element facts.
- **Does correctness depend on it?** No. Proved: the Neumann data `h_j = ½(a_{j−1} − a_{j+1})`
  is the gradient of the area of the boundary polygon with respect to `b`, it sums to zero, and a
  solution of `A b = h` minimizes `½ bᵀAb − 𝒜(a, b)` (`[formula: bff-hilbert-transform]`). The
  quality of the result is measured (`[formula: quasi-conformal-error]`), not assumed. Sign: BFF
  Sec. 4.3.3 solves `A b = −Uᵀa` with `(Uᵀa)_j = ½(a_{j+1} − a_{j−1})` (its Eq. 17 calls this
  `h_j`), so the right-hand side is `½(a_{j−1} − a_{j+1})`, as in the code; its Alg. 6 writes
  `A b = h` with the same `h`, which matches only for the opposite direction of the boundary walk.

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

---

## Out of scope of this package

The end-to-end pipeline (`crates/ddg-uv`: charts, packing), the Python bindings and the Blender
add-on are not covered. For `crates/ddg-bff` and `crates/ddg-optimal-cones` the formulas listed in
the README are covered; their control flow (iteration limits, safeguards, fallbacks, tolerances)
is not.
