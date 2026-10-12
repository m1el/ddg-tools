# Cited, not machine-checked

Results the library relies on, or that the paper uses to motivate it, which are **not** proved in
`lean/proofs`. Each entry has an id (the same `[formula: id]` tag the Rust code may carry;
`scripts/check-formulas.sh` reads the ids from the headings below), an exact reference, the reason
it is not formalized, and whether the correctness of the library *relative to the canonical
discrete problem of `CLAUDE.md`* depends on it.

"Paper" is Soliman, Slepčev, Crane, *Optimal Cone Singularities for Conformal Flattening*,
ACM TOG 37(4), 2018. "Thesis" is Soliman, *Optimal Cone Singularities for Conformal Flattening*,
MS thesis, CMU 2018. Mathlib is the revision pinned in `lake-manifest.json` (tag `v4.34.1`); a
search of that revision finds no subdifferential, no convex conjugate / Fenchel duality, no
Wasserstein distance, no Laplace–Beltrami operator or Gaussian curvature of a Riemannian surface,
and Sobolev spaces only in the form of the Euclidean Sobolev inequality and Bessel potential
spaces.

The short version: **nothing below is needed for the claim "the returned `(u, µ)` solves the
canonical discrete problem"**. That claim rests on the machine-checked sufficiency theorem and
duality-gap certificate. The cited results explain why the canonical discrete problem is the right
problem (smooth theory), why the solver can be expected to terminate (Newton convergence), and
what happens in two optional post-processing steps (cone merging, genus reporting).

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

---

## Out of scope of this package

Boundary First Flattening (Sawhney, Crane 2017), which turns the cone angles into a UV map, and
the cutting of the mesh are implemented in other crates and are not covered here.
