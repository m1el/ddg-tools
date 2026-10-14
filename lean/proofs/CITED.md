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

### [formula: cut-weighted-fem-stiffness] The weighted cotan matrix as a finite element matrix

- **Statement.** For piecewise linear `a`, `f`, `g` on a triangle mesh,
  `∫ a ∇f·∇g = fᵀ L_a g` with `L_a = Σ_t ā_t · (cotan Laplacian of the triangle t)`,
  `ā_t` the mean of `a` over the vertices of `t` (SPEC §3.6, D11).
- **Reference.** Standard finite elements (the gradient is constant on a triangle and the
  integral of a linear function is its mean times the area); the Cuts paper Sec. 5.3 uses the
  edge rule `(a_i + a_j)/2` instead.
- **Why not formalized.** The integral over a triangle is not formalized in this package (see
  `cotan-fem-stiffness`).
- **What depends on it.** Only the name "exact stiffness matrix". Proved without it:
  `xᵀ L_a x = Σ_t ā_t · (cotan energy of x on t) ≥ 0` for `a ≥ 0`
  (`[formula: cut-weighted-laplacian]`), `L_a = cL` for constant `a`, and the edge rule is
  indefinite on an obtuse triangle (`edgeRuleEnergy_neg_example`).

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
is `crates/ddg-intrinsic/SPEC.md` ("SPEC", decisions `N1`–`N18`). "B–S" is A. I. Bobenko,
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
circumcentre in barycentric coordinates. **Not** machine-checked: the six results below; the
halfedge data structure and that the Rust code applies the formulas to the right elements; a
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
- **Why not formalized.** It is a global statement about all geodesic triangulations of a
  surface; neither piecewise flat surfaces nor immersed disks are in Mathlib.
- **What depends on it.** That the result of Delaunay flipping is *the* Delaunay triangulation up
  to ties, independent of the order of the queue (SPEC §5 item 4, N6: tested, the ties are
  counted). Proved without it: the local condition is equivalent to a nonnegative cotan weight
  (`[formula: it-cotan-weight]`, `cotanWeight_nonneg_iff`), and at a tie both diagonals satisfy
  it with equality (`flip_of_tie`).

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
- **Why not formalized.** The argument is about the set of all geodesic triangulations of a
  surface (a function on it that every flip decreases, and finiteness below a bound), which
  rests on the geometry of geodesics on a cone surface.
- **What depends on it.** That `Delaunay flipping` returns at all; the crate has a limit on the
  number of flips and reports when it is reached (SPEC §4, N6; never reached in the recorded
  runs). Proved: every step of the algorithm is possible (an edge that violates the condition has
  a strictly convex diamond, `diamondConvex_of_not_delaunay`; an edge with the same triangle on
  both sides never violates it, `delaunay_of_same_triangle`), and after the flip the new edge
  satisfies the condition strictly (`flip_of_not_delaunay`).

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
- **Why not formalized.** The energy of the interpolant over the two triangulations of a
  quadrilateral could be compared from the cotan formula; it was not asked for and nothing in the
  crate evaluates it. The global statement needs the previous entry.
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
  shorter than the circumradius (`circumradius_gt_of_angle_lt`).

### [formula: it-ruppert-encroachment] Ruppert's boundary rule in the plane

- **Statement.** In the plane, Delaunay refinement that splits a boundary segment whenever a
  vertex lies on or inside its diametral circle and is visible from it (the segment is
  *encroached*), and inserts a circumcentre only if it encroaches on no segment (otherwise the
  circumcentre is rejected and the encroached segments are split), terminates for angle bounds up
  to `arcsin(1/(2√2)) ≈ 20.7°` (circumradius to shortest edge ratio `B ≥ √2`) when no two
  input segments meet at an angle below `60°` (Shewchuk's analysis); it then creates no edge
  shorter than the minimum local feature size of the input. For smaller input angles, splitting
  the segments at an input vertex on concentric circular shells (radii that are powers of
  two) prevents the two segments from encroaching on each other for ever. In a Delaunay
  triangulation, a vertex inside the diametral circle of a boundary segment is opposite an
  obtuse angle (a vertex on it, a right angle), by Thales' theorem.
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
  at `25°` and `30°` is observed, not proved, as in the plane. Proved here: nothing of the rule
  itself; the criterion "the opposite angle is obtuse" is the sign of a coordinate of the
  circumcentre (`circumcenterBary_neg_iff`).

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

---

## Out of scope of this package

The end-to-end pipeline (`crates/ddg-uv`: charts, packing), the Python bindings and the Blender
add-on are not covered. For `crates/ddg-bff` and `crates/ddg-optimal-cones` the formulas listed in
the README are covered; their control flow (iteration limits, safeguards, fallbacks, tolerances)
is not. For `crates/ddg-variational-cuts` the same holds: the formulas of one flow step are
covered, the extraction of the cut and of the cut mesh (which faces are split, the choice of the
quad diagonal, the bookkeeping of regions and components), the run loop around the stop rule,
the `N`-region variant (SPEC §8) as a whole and the helpers for initial level sets are not. For
`crates/ddg-intrinsic` the closed forms of the twelve `it-` ids are covered; the data structure,
the algorithms built on them (flipping, tracing, insertion, overlay, refinement, export) and
their tolerances are not.
