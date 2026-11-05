# Formalization survey of `CITED.md`

Survey (2026-11-05) of the results in `CITED.md` that are cited and not machine-checked: which ones can be proved in
`lean/proofs` with the pinned Mathlib (`v4.34.1`) and the existing `DdgProofs`, which need missing infrastructure,
and which are out of reach. Mathlib names were checked by grep against the pinned sources (and one `#check` file).
Nothing here is proved yet; the classes are estimates.

Classes: **A** provable now (effort S < 300 lines, M 300–1500, L > 1500); **B** needs bounded missing infrastructure;
**C** out of reach (smooth Riemannian/PDE/measure theory beyond Mathlib, or global convergence of a numerical method).
Of the ~116 entries, roughly 30 are A, 27 B and the rest C.

## Findings that concern the code or the claims, not only the proofs

These come first because they may be bugs or wrong statements, independent of whether anyone formalizes anything.

1. **Segment–segment distance in `ccd.rs` can overestimate** (`rc-closest-points`). For nearly parallel lines
   (`0 < denom ≤ 1e-14·a·e`) the interior pair is skipped; an overestimate is the unsafe direction for CCD. Expected
   but unproved bound: `cand² ≤ d² + 1e-14·|q₁−p₁|²`. Prove the bound or fix the code.
   **Status (2026-11-05): confirmed and fixed.** Proved bound `cand² ≤ d² + ¼·1e-14·|q₁−p₁|²`, tight (`segSegThr_sq_le`);
   in f64 the closed form above the threshold was worse (up to ≈3e-6·L near θ = 1e-6, vs 80-digit references), and a
   near-parallel segment could pass through another without a contact. `ccd.rs` now takes `s` from the cross-product
   form and re-projects `t(s)`, `s(t)` (exact in exact arithmetic, f64 error ≈3e-15·L); regression tests in
   `tests/ccd.rs`. The RS CCD uses the same functions.
2. **`rc-accd` termination claim is false in exact arithmetic.** The step `τ += (d−g)/l` can converge to the contact
   time without reaching it; only the `1e-9·d₀` slack makes it finite (up to ~10⁹ steps), so the cap of 100 matters.
3. **`discrete-minimizer-existence` is false as stated**: it fails when the problem is infeasible (closed mesh, every
   vertex forbidden, `χ ≠ 0`). Needs a nonempty-feasible-set hypothesis; whether the Rust code reports infeasibility
   was not checked.
4. **`app-vector-singularity-index` interval mismatch**: CITED says ω ∈ (−π, π), the code uses (−π, π]; an edge at
   exactly π breaks the twin cancellation, so the index sum needs a "no edge at π" hypothesis.
   **Status (2026-11-05): confirmed and fixed.** `Intrinsic/SingularityIndex.lean`: with `arg` on both halfedges
   `Σ_f index = nχ + #{edges with ω = π}` exactly (`sum_faceIndex`); in f64 one exact tie `u_j = −r_h u_i` on each
   halfedge of icosphere(1) gave `index_sum = 3 ≠ 2` in 25 of 240 cases. `singularities` now shifts the higher
   halfedge of a pair by `∓2π` when `|ω_h + ω_t| > π` (tie: `+π`/`−π`; bitwise unchanged elsewhere, pins pass); proved `index_sum = nχ` unconditionally for that rule
   (`sum_round_faceIndexFixed`; closed, oriented, every vertex in a face); `tests/vector_pi_tie.rs` passes (n = 1, 2).
5. **`rc-ncg`** convergence claim needs strong Wolfe conditions; the crate's Armijo-only search does not provide them.

## Ranked candidates (value to library correctness × feasibility)

**Cheap, mostly proved already**

| id | class | what it takes |
|---|---|---|
| `dev-lbfgs-two-loop`, `rc-lbfgs-two-loop` | A-S | already proved as `rs-lbfgs-two-loop` (`RepulsiveSurfaces/Optimize.lean`, `twoLoop_eq`, `bfgsInv_posDef`); instantiate H₀ = γI, and a PSD/subspace variant for RC's mean-free H₀ = γM⁻¹ |
| `app-amr-edge-flux` | A-S | 2 × `divergence_corner` (`IntrinsicApps/HeatFem.lean`) with X = ∇u, plus dot/cross cotangent ↔ `cornerCot`. **Status (2026-11-05): done**, `TriangleIntegral/Flux.lean` (proved directly from `gradInterp` and `cotCorner`, in `side_flux`'s argument order; also twin fluxes add to the jump); moved out of CITED |
| `moreau-yosida-convergence` (finite-dim, with rate) | A-S | `IsRegMinimizer.objective_le` + `IsMinimizer.weighted_sq_sub_le` |
| `discrete-minimizer-existence` (+ feasibility hypothesis) | A-S | compact sublevel set, `IsCompact.exists_isMinOn`; with `IsMinimizer.unique` gives exactly one minimizer |
| `rc-distance-gradient` | A-S | `norm_eq_iInf_iff_real_inner_le_zero`, `HasFDerivAt.sqrt`; closes the `rc-penalty-surface` chain |
| `cut-level-set-evolution` (discrete) | A-S | affine φ on a triangle; discharges the `φ_t = −σ` hypothesis behind the proved cut-flow sign |
| `cut-network-first-variation` (junction part) | A-S | |
| `rc-armijo` (also `dev-armijo-backtracking`), `rc-accd` (with slack) | A-S | differentiability at γ suffices |

**Medium, high value**

| id | class | what it takes |
|---|---|---|
| `rc-closest-points` + `pat-truth-closest-point` | A-M | links the proved `primDist_advance` to the distances `ccd.rs` computes; same lemmas serve the PAT ground truth; settles finding 1. **Status (2026-11-05): done**, `RepulsiveCurves/ClosestPoints.lean` (~1650 lines), both moved out of CITED |
| `semismooth-newton-convergence` (local) | A-M | one Newton step from a nearby iterate lands on the exact zero, plus Stage II nonsingularity; builds on `IsPiecewiseAffine`, `newton_stageII_exact`; the only route to a termination statement |
| `optimality-system-necessary-closed` | A-M | exchange variations + 1-D Helly; with existence gives a KKT triple on closed meshes |
| `dev-simple-eigenvalue-derivative` | A-M | smallest eigenvalue only: proved supergradient + eigenbasis/gap lower bound (`Developability/Covariance.lean`) |
| `dev-min-width-hull-edge`, `dev-max-covariance-voronoi` | A-M / M-L | computed intrinsic max energy equals the minimum, not just an upper bound; the code enumerates all pairs, so no convex hull is needed |
| `rc-near-cos-series` | A-M | Mathlib now has the real-exponent binomial series (`Real.one_add_rpow_hasFPowerSeriesOnBall_zero`) and `hasDerivAt_tsum` |
| `rc-near-closest` (fixed active set) | A-M | the switch points themselves need a Danskin theorem (B-L) |
| `app-steiner-fermat-point` | A-S/M | minimality S; the 120° and S² algebra M; `Triangle.lean`, `Diamond.lean`, `Barycentric.lean` have the trigonometry |

**Shared infrastructure worth building once**

* **Integral over a triangle** (B-M, ~500–1000 lines): affine change of variables to the standard simplex, Fubini, the
  Beta integral for monomials. Covers `cotan-fem-stiffness` (Galerkin vs. lumped mass), `cut-weighted-fem-stiffness`,
  the integral part of `bff-conformal-energy`, `app-vector-fem-connection-laplacian`, `app-fem-galerkin-mass`,
  `it-odt-energy`, and the L² norms of the AMR estimator.
  **Status (2026-11-05): built**, `DdgProofs/TriangleIntegral/` (`Basic.lean` ~450 lines, `Fem.lean` ~710 lines; no
  sorry, axioms the standard three). Per item:
  - `app-fem-galerkin-mass`: **done**, `∫_T λᵢλⱼ = A/6, A/12`, row sums `= ∫λᵢ = A/3`, `∫(Σuλ)² = galerkinForm`.
  - `cotan-fem-stiffness`: **done and moved out of CITED** (`∫∇λᵢ·∇λⱼ = −½cot θ_k`, `∫|∇f|² = cotEnergy`,
    `Σ_f ∫ = uᵀLu` over isometric layouts); the convergence part stays cited as `cotan-fem-convergence`.
  - `app-vector-fem-connection-laplacian`: **partial**, new id `fem-area-coordinate-rule` (rule + the exponential
    slice formula that `monomial_integral` evaluates, tagged in Rust); the gauge model, Eq. 17 and `beta_exp` stay cited.
    **Update (2026-11-05):** `beta_exp` proved as new id `app-vector-beta-exp` (`TriangleIntegral/BetaExp.lean`):
    series with running product and a geometric tail bound, recurrence + binomial expansion, both exact.
    `beta_exp` makes no truncation claim beyond its stopping rule; its rounding remarks are not covered (a numeric
    check of the series one, max_k|t_k| ≤ e^{|Ω|}|B| for p, m ≤ 11, |Ω| < 3, found no violation). Still cited:
    the gauge model and the Eq. 17 closed forms.
  - `cut-weighted-fem-stiffness`: **done, CITED entry removed** (per triangle and `Σ_f ∫ a∇u·∇w = uᵀL_ā w`).
  - `bff-conformal-energy`: **done, CITED entry removed** (2026-11-05): face orientation added
    (`DdgProofs/Orientation.lean`, `TriMesh.Oriented`, twins, boundary halfedges); shoelace over the mesh and
    `2Σ_f∫|∂̄f|² = ½(aᵀLa + bᵀLb) − 𝒜(a, b)` in `Bff/ConformalEnergy.lean`; per triangle `bff-conformal-energy-triangle`.
  - `it-odt-energy`: **partial**, Lemma 4.1 per triangle and for a star as `it-odt-interpolation-error`
    (also `∫_T|u_I − u| = (A/12)Σℓ²`); the tiling of `Ω` and Thm 2.3 stay cited.
    **Update (2026-11-05): done** with `TriangleIntegral/Tiling.lean` (additivity over `IsTiling`, null
    triangle boundaries via `Convex.addHaar_frontier`): Lemma 4.1 over a tiling (also per vertex star and by
    edge lengths), `E(x) = S(x) + C − ∫_Ω|y|²` for a fan plus fixed triangles, and the ODT step lowers `E`
    (tiling of the moved star as a hypothesis). Thm 2.3 and the boundary rule stay cited under the new id
    `it-odt-delaunay-optimality` (needs all triangulations of a point set).
  - AMR estimator: **done for the element norms**, `app-amr-element-norms` (tagged on `amr::estimate`); the bounds
    stay cited. `ic-transfer-mass-integral`: the one-triangle integral is now proved, the tiling stays cited.
    **Update (2026-11-05): `ic-transfer-mass-integral` done and moved out of CITED** (`Tiling.lean`): the
    per-piece assembly of `Transfer::new` is `∫_piece λ_aλ_b`, summing to `∫_T λ_aλ_b` over pieces that tile
    `T` (hypothesis), and `Σ pieceMass = pieceMass_one` for consistently oriented pieces.
  - No disagreement between the cited formulas and the code was found.
* **Multi-face strip trace** (M) for the intrinsic tracing entries.

**Out of reach (C)**: the smooth theory (uniformization, Riemann mapping, fractional Sobolev operators, PDE theory
behind Points as Tori, consistency results of RS/RC, normal isotopy for integer coordinates, global convergence of the
solvers). The "Out of scope" section has nothing in reach.

## Stale or inaccurate reasons in `CITED.md`

* `dev-lbfgs-two-loop`, `rc-lbfgs-two-loop`: already proved (`rs-lbfgs-two-loop`).
* `app-amr-edge-flux`: "cotan-fem-stiffness is cited" is irrelevant; `divergence_corner` and `cotan_formula` cover it. (Resolved: proved, entry removed.)
* `app-steiner-fermat-point`: triangle trigonometry now exists in `DdgProofs`.
* `dev-min-width-hull-edge` ("needs convex hulls"), `dev-max-covariance-voronoi` ("needs the cell decomposition of the
  sphere"): overstated.
* `dev-simple-eigenvalue-derivative`: Mathlib has only the extremal eigenvalue (`hasEigenvalue_iInf/iSup_of_finiteDimensional`),
  no Courant–Fischer.
* `dev-flip-manifold`: `GaussBonnet.lean` models faces, edges, edge degree and χ, enough for χ and degree ≤ 2 under a flip. (Status 2026-11-05: done on gluing data instead, `Delaunay/Flip.lean`; the entry is trimmed.)
* `app-vector-singularity-index`: Gauss–Bonnet is proved (`gauss_bonnet_closed`); the gap is an orientation hypothesis on `TriMesh`.
  **Status (2026-11-05): closed** (orientation added; see finding 4).
* `rc-near-cos-series`: "binomial series with real exponents only in pieces" is stale.
* `rc-distance-gradient`: Mathlib has no projection map onto convex sets, only existence and the characterization.
* `rc-accd`, `rc-ncg`, `discrete-minimizer-existence`: see findings 2, 5, 3.
* `optimal-measure-support`: Jordan and Radon–Nikodym decompositions exist in Mathlib; the obstacle is the smooth optimality system.
* `cone-problem-relaxation`: Riesz–Markov–Kakutani (`rieszMeasure`, positive functionals) and Prokhorov exist; missing are the signed dual and the surface setting.
* `disk-injectivity-criterion`: Mathlib has no winding number at all (nor Jordan curve, Brouwer, degree), not "in pieces".
* `bff-disk-uniformization`: Mathlib now has a partial private `RiemannMapping.lean`; may go stale on the next bump.
* `cone-rounding-stability`: needs no proof; `duality_gap_certificate` already bounds a merged cone configuration after
  the fact, so the Rust side could report that gap.

The per-entry tables follow (one per part of `CITED.md`).


### Survey part 1: CITED.md lines 66–715 (34 entries)

Scope: Smooth formulation (8), Discretization (5), Topology (3), BFF (8), Weights (2), Variational Surface Cutting (8).
Paper abbreviations: OC = Optimal Cones (Soliman–Slepčev–Crane 2018), BFF, MR = Möbius Registration, VSC = Variational Surface Cutting.
"Depends": does library correctness (relative to the canonical problem / SPEC) depend on it. Y / partial / N.
Classes: A = provable now (S <300, M 300–1500, L >1500 lines); B = needs bounded missing infra; C = out of reach.
Mathlib paths are relative to `lean/proofs/.lake/packages/mathlib/Mathlib`. I checked them by grep. I compiled nothing.

#### Table

| # | id | paper | depends | class / effort | partial option | key Mathlib / DdgProofs decls |
|---|----|-------|---------|----------------|----------------|-------------------------------|
| 1 | smooth-yamabe-equation | OC | N | C | Flat case only: `log‖f'‖` is harmonic for holomorphic `f` (K=0). Low value | Mathlib: `AnalyticAt.harmonicAt_log_norm` (Analysis/InnerProductSpace/Harmonic/Constructions.lean:109), `InnerProductSpace.Laplacian` (Euclidean only). No curvature. Discrete analogue already proved |
| 2 | cone-problem-relaxation | OC | N | C | The finite-dim analogue is #12 (existence) | RMK exists for positive functionals (`rieszMeasure`, `integral_rieszMeasure`, MeasureTheory/Integral/RieszMarkovKakutani/Real.lean). Prokhorov exists. No signed `C(M)* ≅ M(M)`, no surfaces |
| 3 | smooth-predual-duality | OC | N | C | Already done discretely (optimality-system-sufficient/-necessary, duality-gap-certificate) | 0 hits for subgradient/subdifferential/convexConjugate/Fenchel. DdgProofs has its own `IsSubgradientWithinAt` (SoftThreshold.lean:33) and `Site.conj` (Duality.lean) |
| 4 | poisson-measure-data | OC | N | C | Discrete PD/PSD already proved | Sobolev = Bessel potential only: `TemperedDistribution.MemSobolev` (Analysis/Distribution/Sobolev.lean) and Analysis/FunctionalSpaces/SobolevInequality.lean. No elliptic regularity |
| 5 | optimal-measure-support | OC | N | C | Discrete version done (`eq_softThreshold_add_iff_cases`) | Jordan decomposition and RN do exist: `SignedMeasure.toJordanDecomposition` (MeasureTheory/VectorMeasure/Decomposition/Jordan.lean:200), …/RadonNikodym.lean. See stale list |
| 6 | optimal-measure-regularity | OC | N | C | none (discretely false/meaningless: clusters occur) | — |
| 7 | cone-rounding-stability | OC | partial (optional merge only) | C (no Wasserstein: 0 hits) | **No proof needed, a runtime certificate instead:** the merged pair is exactly feasible, and the pre-merge dual point is still feasible. So `duality_gap_certificate` (Duality.lean:233) / `adjoint_gap_certificate` (Duality.lean:289) bounds the merged suboptimality a posteriori. Rust could report `J(merged) − D` | `duality_gap_certificate`, `weak_duality` |
| 8 | moreau-yosida-convergence | OC | N | **A, S (<100)** | Finite-dim with a rate: `¼ Σ wE m (u_γ − u*)² ≤ Σ γ_i/2 µ*_i²`, hence `µ_γ − µ* = L(u_γ − u*)` → 0 (O(√γ)) when all `wE m > 0` | `IsRegMinimizer` (Newton.lean:394, includes `Feasible`), `IsRegMinimizer.objective_le` (Newton.lean:566), `IsMinimizer.weighted_sq_sub_le` (Optimality.lean:312) |
| 9 | variational-discretization | OC | N | C (FE over measure spaces) | none worth it | — |
| 10 | cotan-fem-stiffness | OC | N (naming) | B, M (500–1000) | Build **triangle-integration infra**: ∫ over a 2-simplex of affine/polynomial functions. That gives (a) `∫_T|∇f|² = Area·|∇f|²`, (b) the Galerkin mass `Area(1+δ_ij)/12` and its row sum = the lumped mass `Area/3`. Convergence stays C | Mathlib has no simplex volume (no `volume (convexHull …)`/stdSimplex hits). Use `volume_regionBetween_eq_integral`, Fubini, `integral_image_eq_integral_abs_det_fderiv_smul` (MeasureTheory/Function/Jacobian.lean:1213), `addHaar_image_continuousLinearMap` (Measure/Lebesgue/EqHaar.lean:310). DdgProofs: `cotan_formula` (Cotan.lean:118) |
| 11 | semismooth-newton-convergence | N (non-termination detected at run time), but it is the "not covered at all" item | global: C (can cycle; mesh independence is ∞-dim). **Local: A, M (400–800)** | (i) **Local exact one-step termination:** if `x*` is an exact zero, the Newton matrix at `x` is injective, and each `z_i = φ_i+µ_i` lies in a piece whose closure contains `z*_i` (radius = distance from `z*_i` to the nearest breakpoint ≠ `z*_i`), then one Newton step lands exactly on `x*`. Proof: `F(x)+DF(x)(x*−x) = F(x*) = 0`, then uniqueness. (ii) Stage II Newton matrix nonsingular for symmetric invertible `L` and `wE m > 0` (complementary supports of `ν`, `q` ⇒ `⟨Lv,q⟩ = −Σ w v² = 0`) | DdgProofs `IsPiecewiseAffine`, `newton_stageII_exact`, `residual_eq_zero_iff`, `residualDeriv` (Newton.lean). Nonsingularity not proved anywhere yet. Mathlib: `Matrix.PosDef`, `Matrix.IsHermitian` |
| 12 | discrete-minimizer-existence | N (but upgrades "unique" to "∃!"; with `existsUnique_kkt` it gives an existing KKT triple) | **A, S (150–250)** if all `wE_i m_i > 0` and feasible; general case B, L | Statement must assume feasibility (see stale list). Proof: minimize `u ↦ J(u, Lu+Ω)` on the closed set `{u : (Lu+Ω)_i ∈ feasible_i}`; the sublevel set is bounded, so compact. Add feasibility lemmas: Dirichlet always feasible (`µ` = 0 / fixed angles); closed connected feasible ⟺ ∃ admissible `µ` with `Σµ = ΣΩ`. General case needs closedness of the linear image of a polyhedron (Minkowski–Weyl / Fourier–Motzkin) | `IsCompact.exists_isMinOn`, `ContinuousOn.exists_isMinOn'` (Topology/Order/Compact.lean:252), `Metric.isCompact_of_isClosed_isBounded` (Topology/MetricSpace/Bounded.lean:320), `continuous_coe_real_ereal` (Topology/Instances/EReal/Lemmas.lean:52). DdgProofs: `IsMinimizer.unique`, `IsMinimizer.existsUnique_kkt` (Optimality.lean:324, 351), `TriMesh.exists_laplacian_mulVec_eq_iff` (Bff/Laplace.lean:185). Minkowski–Weyl is only *mentioned* in a docstring (Geometry/Convex/Cone/DualFinite.lean:21); there is no proof |
| 13 | optimality-system-necessary-closed | N (converse only) | **A, M (400–900)** | Pairwise exchange `µ_j += t, µ_k −= t` keeps `Σµ`. One-sided derivative bounds give `a_j + y_j ≤ b_k + y_k` for all `j,k`. A 1-D Helly step gives a common `ν`. Then `φ` from range `L = 1^⊥`. Also closes the FactorOnce closed converse | DdgProofs `stationary_of_minimizer` (Optimality.lean:79), `IsMinimizer.weighted_mean_eq_zero` (Closed.lean:111), `exists_laplacian_mulVec_eq_iff`, FactorOnce `pinBlock`, `reducedOptimal_of_isMinimizer_dirichlet` (FactorOnce.lean:440). Mathlib `LinearMap.IsSymmetric.orthogonal_range` (Analysis/InnerProductSpace/Symmetric.lean:332) |
| 14 | genus-from-euler | partial (meaning of the `is_disk` check) | C | `χ ≤ 2C − b` (g ≥ 0) by tree–cotree counting: B, M–L (needs #15(i)). Integrality of g needs orientability plus homology: C/L+ | No classification of surfaces. Mathlib has only `AlgebraicTopology/SimplicialComplex`, `Algebra/Homology/EulerCharacteristic.lean`. DdgProofs `TriMesh.eulerChar`, `Manifold` (GaussBonnet.lean:101, 105) |
| 15 | tree-cotree-disk | N (charts checked at run time) | (i) B, M (600–1200); (ii) C (disk topology) | (i) is a pure parity argument, see notes. Combinatorial part of (ii) (cut mesh has χ=1, one boundary loop): M–L, but it only duplicates the runtime check | `SimpleGraph.Connected.exists_isTree_le`, `IsTree.exists_vert_degree_one_of_nontrivial`, `isAcyclic_iff_forall_isBridge`, `IsTree.card_edgeFinset` (Combinatorics/SimpleGraph/Acyclic.lean). DdgProofs `faceEdges`, `edgeDegree`, `boundaryEdges` (GaussBonnet.lean) |
| 16 | cut-leaf-pruning | N | χ part done; loop count/manifoldness B, M–L; disk C | Needs cyclic boundary-loop infrastructure (absent in DdgProofs). Low value | `cut_euler_characteristic` (Bff/Cut.lean:311) |
| 17 | bff-target-boundary-length | defines the target | n/a (a definition) | Possibly a one-liner tying it to the edge lengths of `TriMesh.conformalScale`; otherwise nothing | Yamabe.lean, Bff/Curve.lean |
| 18 | bff-smooth-cherrier-formula | N | C on surfaces | Planar case: for holomorphic `f`, `f' ≠ 0`, `u = log|f'|`: `Δu = 0` and `κ̃ = e^{−u}(κ + ∂u/∂n)` along a parametrized curve. B, M. Motivation only | `AnalyticAt.harmonicAt_log_norm` |
| 19 | bff-conformal-energy | N (quality measured) | **A, S–M (200–600)** | Per triangle: `Area·|∂̄f|² = ¼·Area(|∇a|²+|∇b|²) − ½·sda(f(T))`, and it is 0 iff `f|_T` is an orientation-preserving similarity. Summing the image areas to the loop area `𝒜` needs a consistent-orientation predicate on `TriMesh` (absent), +M. `∫` identification uses #10's infra | `cotan_formula` (Cotan.lean:118), `sda` (Bff/Quality.lean:132), `loopSignedArea` (Bff/Curve.lean:989), `hilbertTransform_minimizes_conformalEnergy_iff` |
| 20 | bff-disk-uniformization | N | C | none (the iteration is a heuristic) | Analysis/Complex/RiemannMapping.lean has only *private* partial steps; the main theorem is not merged (PR #33505) |
| 21 | disk-injectivity-criterion | partial (meaning of `ChartReport::is_injective`) | B, L (2500–4000) | **Discrete degree lemma** (L, ~1500): on a consistently oriented mesh with all images positively oriented, `#{T : p ∈ f(T)}` = signed ray-crossing count of the boundary polygon for generic `p` (interior edges cancel). The polygonal Jordan step "simple ⇒ winding ≤ 1" is more L | Mathlib: 0 hits for winding number, Jordan curve, Brouwer FPT, mapping degree. DdgProofs `sda`, segment-proper-crossing |
| 22 | mobius-centering-existence | partial (canonical up to rotation) | C | none cheap | no Brouwer/degree in Mathlib. Done: `sphereInversion_mem_sphere`, `hasDerivAt_sphereInversion` (Bff/Sphere.lean) |
| 23 | spherical-triangle-area | partial (`covering_degree`, sphere flag) | area meaning (Girard) C; **algebra A, M (300–800)** | (a) `tan(E/2) = det/(1+a·b+b·c+c·a)` with `E = A+B+C−π`, `A,B,C` from the dihedral-angle formula; or (b) cocycle additivity mod 2π (what a degree sum needs), via `Quaternion` | `Measure.toSphere` (MeasureTheory/Constructions/HaarToSphere.lean:51) exists but there are no geodesic-polygon areas. `InnerProductGeometry.angle` |
| 24 | sphere-injectivity-criterion | partial | C | Sphere analogue of #21's degree lemma, L, after #23(b) | — |
| 25 | feature-size-weights | N | A, S (<30) | `r_i > 0` under both readings for `ε > 0`, `m > 0` (the CITED claim). Nothing else to prove | cf. `cut-feature-size` (Cuts/Step.lean:678) |
| 26 | cone-feedback-rule | N | A, S–M (200–400) | Abstract state machine: weights nonincreasing; each vertex reweighted at most once; at most `max_rounds` solves; every solve is `WellFormed` with a positive weight per closed component. Today these hold "by reading" | ConeFeedback.lean (mismatch only), `ConeProblem.WellFormed` |
| 27 | cut-yamabe-equation | N (modelling) | C | Done on closed forms (`diskU_laplacian`, `capU_ode`) | — |
| 28 | cut-shape-derivative | partial (flow direction) | C on surfaces; planar B, L | Star-shaped planar domains in polar coordinates: M–L. Poor value/cost | divergence theorem only on boxes: `integral_divergence_of_hasFDerivAt_off_countable` (MeasureTheory/Integral/DivergenceTheorem.lean:266); `integral_comp_polarCoord_symm` (Analysis/SpecialFunctions/PolarCoord.lean:138) |
| 29 | cut-cea-method | partial (flow direction) | C | Finite-dim Céa (adjoint = derivative w.r.t. a parameter of `A(θ)u = f(θ)`): A, S. Buys little, because the code's density is not the derivative of the discrete objective in the cut position (the combinatorics change) | `hasFDerivAt_ringInverse` (Analysis/Calculus/FDeriv/Mul.lean:661); DdgProofs `scaledEnergy_perturb` |
| 30 | cut-level-set-evolution | partial: the flow's sign rests on `φ_t = −σ` | Euclidean: **A, S** (chain rule) / M (curve existence via IFT); surfaces C | **(a) Discrete:** for `φ` affine on a triangle (as the code interpolates), the zero line moves with normal speed `−φ_t/|∇φ|`; that is `−φ_t` for a signed distance. S (~150). (b) Smooth planar: a curve on the zero set has `φ_t = −∇φ·γ'`; existence of the curve by the IFT. Together they discharge the hypothesis of `firstVariation_levelSet` | `implicitFunctionOfBivariate` (Analysis/Calculus/ImplicitFunction/Bivariate.lean:49), `HasFDerivAt.comp`; DdgProofs `crossingParam`, `interp_crossingParam` (Cuts/Crossing.lean) |
| 31 | cut-ripple-linearization | partial (length-weight choice) | linearization C; missing ingredient A, M (300–600) | Curvature `(m²−1)ε/R²` of the rippled equator, to first order, via `κ_g = det(γ,γ',γ'')/|γ'|³` on S² | DdgProofs Cuts/Ripple.lean, `sphere_mode_ode` |
| 32 | cut-network-first-variation | partial (2π/3 rule, non-default stages) | junction part **A, S**; smooth curve part B, M | (a) `∇_p Σ w‖p−q_k‖ = w Σ unit(p−q_k)`; three unit vectors sum to 0 ⟺ pairwise angles 120°. (b) The stated hypothesis "angles at a triple point sum to 2π" via oriented angles. (c) First variation of length `= −∫κσ` by integration by parts | `contDiffAt_norm`, `hasStrictFDerivAt_norm_sq` (Analysis/InnerProductSpace/Calculus.lean:154, 197), `Orientation.oangle_add_cyc3` (Geometry/Euclidean/Angle/Oriented/Basic.lean:474), `intervalIntegral.integral_mul_deriv_eq_deriv_mul`. DdgProofs `junction_loads_sum` (Cuts/Regions.lean) |
| 33 | cut-weighted-fem-stiffness | N (naming) | B; **S once #10's infra exists** | `∫_T a∇f·∇g = Area·ā·∇f·∇g` for linear `a` | DdgProofs `weightedLaplacian`, `faceMean` (Cuts/Weighted.lean) |
| 34 | cut-fast-marching-convergence | N | convergence C (viscosity solutions); causality A, M (300–500) | Acute-triangle causality: with the minimizing λ in (0,1) and a non-obtuse angle at `i`, `T_i ≥ max(T_j,T_k)`. Single-pass = fixed point: L (needs an algorithm model). Low value | DdgProofs `hopfLaxTriangle_isLeast`, `min_le_hopfLaxTriangle` (Cuts/Redistance.lean) |

#### Notes (non-obvious entries only)

- **#12 existence.** Use the variable `u` only. `µ := Lu + Ω` is determined by `u`, and the feasible set `S = {u | ∀i, (Lu+Ω)_i ∈ (site i).feasible}` is closed: `boundSet` is a preimage of closed sets under `continuous_coe_real_ereal`, and `fixed a` gives a singleton.
  - On `S`, `J ≥ ½ c‖u‖²` with `c = min wE m > 0`, so `{u ∈ S | J ≤ J(u₀)}` is compact. The proof works for **every** `L`, closed or Dirichlet, so it is uniform with `IsMinimizer.unique`.
  - Result: `Feasible nonempty ∧ ∀i, 0 < wE_i m_i ⇒ ∃! minimizer`. With `existsUnique_kkt`, the Dirichlet KKT system has a solution. With #13, so does the closed one.
- **#13 closed necessity.** The "common point of per-vertex subdifferential intervals" is a 1-D Helly argument: pairwise `max_j(a_j+y_j) ≤ min_k(b_k+y_k)`, where `a_j, b_j` are the one-sided slopes of `λ·cost_j + ι_feasible` at `µ_j`. These intervals are already characterized in SoftThreshold.lean (`isSubgradientWithinAt_abs_iff`).
- **#11.** No Mathlib semismooth calculus is needed. The repo's `IsPiecewiseAffine` is exactly the right abstraction. The new hypothesis is "piece compatibility", `Q(z*) = Q(z) + slope(z)(z* − z)`, which follows from continuity of `Q` and `z*` lying in the closure of `piece(z)`.
  - For soft-threshold and clamp sites this is a concrete radius: `|z − z*| <` the distance to the nearest breakpoint `≠ z*`.
  - This turns "termination not covered at all" into a proved local statement: any iterate within that radius terminates in one step, given the nonsingularity lemma.
- **#15(i).** Let `S` be a face set closed under dual adjacency outside `T*`. At every vertex `v`, double counting face/edge incidences in the star gives: `#(∂S edges at v) + #(boundary edges at v whose face is in S)` is even.
  - Contracting the boundary loops therefore makes `∂S ⊆ T` a nonempty even edge set, which contradicts acyclicity.
  - No vertex-link connectivity is needed: the existing `TriMesh.Manifold` plus "every interior edge in exactly 2 faces" suffices.
  - Define acyclicity of `T` in cycle-space form ("no nonempty even subset"). That avoids multigraphs and loops from the contraction; Mathlib's `SimpleGraph` can't represent those.
- **#14/#15/#16, what the topology buys.** The BFF solves need only connectivity, one boundary loop and Gauss–Bonnet, all of which are proved. Genuine disk-ness matters only for the meaning of the injectivity flag (#21). Priority for this block is low.
- **#10/#33/#19 (shared infra).** One file, "integrals over a triangle in ℝ²", gives:
  - `∫_T c = Area`;
  - `∫_T` affine `=` mean·Area;
  - `∫_T λ_iλ_j = Area(1+δ_ij)/12`.

  Route: Fubini on the reference triangle `{s,t ≥ 0, s+t ≤ 1}` (`regionBetween`), then an affine change of variables. It is also reused by the other share's FEM entries (CITED ~l.1176–1191: the barycentric monomial rule `∫_T b^α = 2|T| α!/(|α|+2)!`). It is cross-cutting, so its value is higher than any single entry suggests.
- **#7.** Prefer a runtime certificate over Mainini's `H^{-2}`/W₂ estimate. A merged configuration is exactly feasible whenever the merged site admits the summed angle, so the existing gap theorems already bound how much worse it is.

#### Ranked top candidates (value to library correctness × feasibility)

1. **discrete-minimizer-existence, case `wE m > 0` + explicit feasibility** (A, S). Fills the "existence not proved" hole that CITED's preamble names, upgrades uniqueness to ∃!, and yields existence of the KKT triple (Dirichlet) the Newton solver targets.
2. **semismooth-newton local exact termination + Stage II Newton-matrix nonsingularity** (A, M). The only route to any proved statement about termination, the item the preamble lists as "not covered at all". It builds directly on Newton.lean.
3. **moreau-yosida-convergence, finite-dim with rate** (A, S, ~50–100 lines). Two existing lemmas combine to `¼Σ wE m (u_γ−u*)² ≤ Σ γ_i/2 µ*_i²`, which quantifies the Stage I warm start.
4. **optimality-system-necessary-closed** (A, M). Completes the KKT ⟺ minimizer iff on closed meshes. Combined with #1, it gives existence of a KKT triple on closed meshes, and it also closes the FactorOnce closed converse.
5. **Triangle-integration infrastructure → cotan-fem-stiffness (a,b), cut-weighted-fem-stiffness, bff-conformal-energy `∫` identification** (B, M). Cross-cutting, and also serves the FEM entries of the other shares.
6. **cut-level-set-evolution, discrete affine version (+ planar chain rule)** (A, S). Discharges the hypothesis `φ_t = −σ` on which the proved sign of the cut flow (`firstVariation_levelSet`) rests.
7. **bff-conformal-energy per-triangle LSCM identity** (A, S–M).
8. **cut-network-first-variation junction part + the 2π angle-sum hypothesis** (A, S). It removes a stated unproved hypothesis of `junction_loads_*`.
9. **tree-cotree-disk (i)** (B, M). Pure parity plus Mathlib spanning trees.
10. **disk-injectivity discrete degree lemma** (B, L). Gives the injectivity flag a partial meaning and is reusable for the sphere.

Also cheap but of near-zero value: feature-size-weights positivity (S), and cone-feedback-rule invariants (S–M, would replace "by reading").

#### Stale, incorrect or imprecise CITED.md claims

1. **discrete-minimizer-existence (l.213–221): the statement is false as written.** "The canonical discrete problem has a minimizer" fails when the problem is infeasible. Example: a closed mesh with every vertex forbidden (`µ = 0`) and `χ ≠ 0`, since `Σµ = ΣΩ = 2πχ` is then impossible. Bounds that are too tight do the same. The objective is always bounded below by 0, so "bounded below on a polyhedron" hides the real hypothesis, a nonempty feasible set. Restate as "if feasible, a minimizer exists". (Whether the Rust side reports infeasibility explicitly is worth checking; that is outside my share.)
2. **optimal-measure-support (l.127–128):** the reason names "Jordan decomposition, Radon–Nikodym" as the missing tools, but both are in Mathlib (`MeasureTheory.SignedMeasure.toJordanDecomposition`, MeasureTheory/VectorMeasure/Decomposition/{Jordan,RadonNikodym}.lean). The real obstacle is the smooth optimality System 11 on a surface (i.e. smooth-predual-duality).
3. **cone-problem-relaxation (l.87):** "Needs `M(M) = C(M)*` (Riesz–Markov)". Mathlib *has* Riesz–Markov–Kakutani for positive functionals on `C_c(X,ℝ)` (`rieszMeasure`, `integral_rieszMeasure`), and Prokhorov. What is missing is the signed-dual identification plus the surface setting. Reword.
4. **disk-injectivity-criterion (l.394–395):** "winding numbers of polygons… which Mathlib has only in pieces" overstates it. The pinned Mathlib has **no** winding number, Jordan curve theorem, Brouwer fixed point or mapping degree (grep: 0 hits).
5. **semismooth-newton-convergence (l.205):** "Mathlib has no Newton (semismooth) differentiability" is true but not the obstacle. The repo's own `IsPiecewiseAffine` suffices for local finite termination (#11). "Not attempted" is the accurate reason.
6. **moreau-yosida-convergence (l.162–163):** "the finite-dimensional analogue was not attempted" is accurate, but it is a few lines from `IsRegMinimizer.objective_le` + `IsMinimizer.weighted_sq_sub_le`, so it should not stay cited.
7. **cone-feedback-rule (l.494–498):** "nothing to prove about a threshold" is true of the threshold. The termination and well-formedness claims ("by reading") are formalizable invariants.
8. **bff-disk-uniformization (l.375):** "Uniformization is not in Mathlib" is still true at the pin. However, Analysis/Complex/RiemannMapping.lean now holds private partial steps toward the Riemann mapping theorem (upstream PR #33505), so this may go stale on a Mathlib bump. It would not help with the iteration anyway.
9. Header (l.14–18) is accurate for the pin: no subdifferential, conjugate, Fenchel, Wasserstein, Laplace–Beltrami or Gaussian curvature, no classification of surfaces, and Sobolev only as Bessel potential spaces plus the Sobolev inequality. It could add that RMK, Jordan/RN decomposition, the implicit function theorem (`implicitFunctionOfBivariate`) and SimpleGraph spanning trees are available.

### CITED.md survey, part 2 (lines 716–1494): Intrinsic, Apps 5.1–5.5, Developability

Pinned Mathlib: `lean/proofs/.lake/packages/mathlib/Mathlib` (v4.34.1). Paths below are relative to that
directory (Mathlib) or to `lean/proofs/DdgProofs` (repo). Nothing was compiled; everything was found by grep.
29 entries. Depends: Y = a formula or a claim the library makes rests on it; partial = only the
"min, not just an upper bound" or "terminates" part; N = motivation or reporting only.

#### Table

| id | paper | depends | class / effort | partial option | key Mathlib / repo decls |
|---|---|---|---|---|---|
| it-delaunay-existence-uniqueness | B–S 2007 | partial (order-independence of flips, N6; tested) | C | none useful beyond what is already proved (`cotanWeight_nonneg_iff`, `flip_of_tie`) | no Delaunay, no piecewise flat surfaces in Mathlib (grep: 0 hits) ; **Status (2026-11-05): layer 1** — B–S Lemma 9 for one diamond (empty circumcircle ⟺ θk+θl ≤ π) proved, `Delaunay/Local.lean`; the global statement stays cited |
| it-delaunay-flip-termination | Indermitte+ 2001, B–S | partial (the flip loop is capped, N6) | C (surface); planar case B-L (done in Coq by Dufourd–Bertot) | **B-M**: one-diamond strict decrease of a potential (Musin harmonic index Σ(a²+b²+c²)/A, or Σ circumradius²) when θk+θl>π, as a 5-length inequality | repo `Intrinsic/Diamond.lean` (`flippedLengthFromLengths_*`, `flip_isTriangle`, `flip_area`), `Intrinsic/Delaunay.lean` (`flip_of_not_delaunay`), `Barycentric.lean` (`circumradius_eq_div_sin`) ; **Status (2026-11-05): one-flip part done**, `Delaunay/Potential.lean` (`it-delaunay-potential`): Σcircumradius² was not used; the lifting energy ΣA(a²+b²+c²) changes by exactly −8A_ijk A_jil(cot θk+cot θl), and the harmonic index follows from Rippa; termination (finiteness) stays cited; **layer 2 (2026-11-05)**: the flip on a whole triangulation (gluing data with self-loops etc., `Delaunay/{Complex,Flip,Metric}.lean`): every non-Delaunay interior edge between two distinct faces can be flipped, the result is again a metric triangulation with the same V, E, F, χ, boundary and cone angles; **layer 3 (2026-11-05)**: `Delaunay/{Termination,FlipQueue}.lean` (`it-delaunay-self-glued`, `it-delaunay-flip-algorithm`, `it-delaunay-flip-loop`): self-glued edges are Delaunay, global potentials decrease, termination + Delaunay end state + invariants under an explicit finiteness hypothesis (finitely many triangulations with bounded edges up to isomorphism, the only part still cited), Gauss–Bonnet with boundary, and the Rust queue loop (tolerance, any order) |
| it-rippa-theorem | Rippa 1990, B–S | N (N6 remark, tested) | A-M for one diamond; global part C (needs the previous entry) | E_T1(f)−E_T2(f) = κ·ℓ(f)², κ ≤ 0 iff the new edge is Delaunay (rank ≤ 1 because both energies agree on affine f and share the boundary term) | repo `Cotan.lean` (`cotan_formula`, `cotEnergy_identity`), `Diamond.lean` layouts ; **Status (2026-11-05): one diamond done**, `Delaunay/Rippa.lean` (`it-rippa-diamond`), κ = −(cot θk+cot θl)/(2A_ilk A_jkl), ℓ = affine dependence with the four triangle areas; proved as a coordinate identity, not by the rank argument; global part cited |
| it-chew-second-algorithm | Chew 1993, Shewchuk | N (N1 says "argued/observed") | C (planar Delaunay refinement, packing argument) | ingredients already proved (`circumradius_gt_of_angle_lt`) | none in Mathlib ; **Status (2026-11-05)**: R = a/(2 sin θ), θ < π/6 ⇒ a < R also from lengths only (`Delaunay/Encroachment.lean`) |
| it-ruppert-encroachment | Ruppert 1995, Shewchuk | N | C | A-S: "inside the diametral disk iff the angle is obtuse" (Thales, strict form) | `Geometry/Euclidean/Angle/Sphere.lean` (`thales_theorem`), `inner_neg_iff_pi_div_two_lt_angle` (`Angle/Unoriented/Basic.lean`); repo `circumcenterBary_neg_iff` ; **Status (2026-11-05): A-S done**, `Delaunay/Encroachment.lean` (`it-encroachment-angle`, tagged on `refine.rs` `is_encroached`); termination cited |
| it-straightest-geodesic | Polthier–Schmies 1998 | Y (a trace evaluates exp; all signpost queries) | B-M | **multi-step trace = one straight line in the strip unfolding** (induction over a list of crossed faces given by lengths, composing the proved one-step lemmas with rigid motions); vertex case stays cited (N7 stops there) | repo `Intrinsic/Trace.lean` (`exit_cases`, `exit_ahead_*`), `Diamond.lean` (`unfoldApex_eq_rigid`, `dist2_rigid`, `sda_rigid`), `IntegerCoords/Strip.lean` (`hit_point`, `hit_unique`) |
| it-odt-energy | Chen–Holst 2011 | partial (N21; the planar star descent is proved) | B-M for Lemma 4.1 (S = ∫|u_I−u|); Thm 2.3 C | Lemma 4.1 on one triangle needs "∫ over a triangle" infra (see shared infra below) | `MeasureTheory/Function/Jacobian.lean` (`integral_image_eq_integral_abs_det_fderiv_smul`), `MeasureTheory/Integral/Prod.lean` (`setIntegral_prod`), `integral_pow` |
| it-odt-target-length-splits | Tournois+ 2009 | N (heuristic schedule) | — (no statement) | none | — |
| app-steiner-fermat-point | Gilbert–Pollak, Kimberling | Y (place and length of every Steiner point, N30) | A-M | (i) **A-S**: if the three unit vectors from F to the corners sum to 0 (⇔ 120°), F minimizes Σ‖x−pᵢ‖ (one Cauchy–Schwarz line); (ii) A-M: the barycentrics aₖ/sin(θₖ+π/3) (∝ 1/(4T+√3(b²+c²−a²))) give dⱼ²+dₖ²+dⱼdₖ = aᵢ² and Σd = S with the S² formula; network minimality B-M | `law_cos` (`Geometry/Euclidean/Triangle.lean`), `EuclideanGeometry.angle`; repo `Triangle.lean` (`cos_cornerAngle`, `sin_cornerAngle`, `angle_eq_cornerAngle`), `Barycentric.lean` (`cornerDistances_eq_baryLength`) |
| app-vector-transport | VHM Eq. 10, GODF Eq. 12 | Y, but it is a definition | — (definition) | could state "ρ is the rotation that maps edge direction at i to edge direction at j in the diamond layout" (A-S), low value | repo `Intrinsic/AppsAmrVector.lean` |
| app-vector-fem-connection-laplacian | GODF Sec. 6, App. D | Y for `ConnectionLaplacian::FiniteElement` only | B-M/L | **area-coordinate rule** ∫_T b₀ᵃb₁ᵇb₂ᶜ = 2|T|a!b!c!/(a+b+c+2)! (what Rust uses) and the Eq. 17 closed form, given the integrand; the gauge model itself stays a modelling choice | `Gamma/Beta.lean` (`betaIntegral_eval_nat_add_one_right`), Jacobian change of variables, Fubini as above |
| app-vector-singularity-index | GODF Eq. 19, App. B | N (reported index sums) | per-face integrality A-S; Σ = nχ B-M | per face: from `e^{iω}r` products and `transport_around_face`; global: needs an *oriented* closed-mesh hypothesis (each directed edge in ≤1 face) on top of `TriMesh.Manifold`, then `gauss_bonnet_closed` | `Complex.exp_eq_one_iff`; repo `GaussBonnet.lean` (`gauss_bonnet_closed`, `eulerChar`), `AppsAmrVector.lean` (`transport_around_face`, `scaled_minus_sum`) |
| app-vector-maximum-principle | VHM Sec. 5.4 | N (a check routine) | A-S (~30 lines) | u_i = Σ wᵢⱼ rⱼᵢ uⱼ / Σw with w ≥ 0 is a nonnegative combination | repo `cotanWeight_nonneg_iff` |
| app-amr-edge-flux | cotan formula | Y (jump residuals of `amr::estimate`) | **A-S** (~100–200 lines) | full statement: it is 2·`divergence_corner` with X = ∇u (`gradSum_linear`), plus cot(dot/cross) = `cornerCot` from lengths | repo `IntrinsicApps/HeatFem.lean` (`divergence_corner`, `gradSum_linear`), `Triangle.lean` (`cornerCot_eq_cot`), `Cotan.lean` (`cot_mul_area`, private) |
| app-amr-residual-estimator | MN 2005, MNS 2002 | N (marking only; no claim about solves) | C (H¹ on polygons, Clément interpolation) | the per-face L² norms of linear R as Galerkin-mass quadratic forms (A-S, after the ∫-over-triangle infra) | Mathlib Sobolev only via Fourier/Bessel potentials (`Analysis/Distribution/Sobolev.lean`, `BesselPotentialSpace.lean`); no domain H¹₀ |
| app-amr-bulk-marking | MNS Sec. 3, MN | N | contraction C; greedy A-S | the greedy loop of `amr::mark` returns a set carrying the fraction θ² (≥, or > for `Indicators`) of minimal cardinality | `Finset.sum`, `List.mergeSort` lemmas |
| dev-hinge-coplanar | Dev. Thm 3.3 | N (interpretation of energies) | A-M with algebraic stand-ins for embeddedness | hinge ⇒ coplanar is trivial; coplanar ⇒ hinge needs only "no two edge vectors of the star positively parallel" and "no fold (adjacent normals not antiparallel)", then pigeonhole on cyclic changes of normal | repo `Developability/Vec.lean` (`cross_cross`, `dot_cross_self_*`, `norm2_cross_add_sq_dot`) |
| dev-hinge-antiparallel | Dev. Prop A.1 | N | A-M, same hypotheses | as above | as above |
| dev-hinge-ruled | Dev. Prop A.2 | N | B-L (global path argument on a mesh) | none worth it | — |
| dev-valence-three-flat | Dev. Prop A.3 | N (default `omit_valence_three`) | A-S given A.1 | two antiparallel edges of one triangle ⇒ zero normal | repo `Vec.lean` |
| dev-simple-eigenvalue-derivative | Magnus 1985, Kato | **Y** (gradient of E_λ and intrinsic width at a simple λ_min) | **A-M** (~400–700 lines) | for the *smallest* eigenvalue no perturbation theory is needed: upper bound = proved supergradient; lower bound from the eigenbasis expansion and the gap g: λ_min(A+H) ≥ λ + x₀ᵀHx₀ − O(‖H‖²/g); squeeze gives `HasDerivAt` of λ_min(A(t)) | repo `Developability/Covariance.lean` (`IsEigenbasis`, `quad_eq_of_eigenbasis`, `quad_isLeast`, `eigenvalue_supergradient`, `hasDerivAt_quad_covMatrix`); Mathlib `Rayleigh.lean` (`hasEigenvalue_iInf_of_finiteDimensional`), `Matrix.IsHermitian.eigenvectorBasis` |
| dev-max-covariance-voronoi | Dev. App. B.4 | partial (λ_max is the min, not only an upper bound) | A/B-M-L (~800–1500) | no spherical Voronoi needed: at a global min (exists, compact sphere) the active signed sites number ≥ 3 distinct points (1 or 2 active ⇒ a linear functional on a sphere or a great circle has no positive local min), 3 distinct unit vectors are never collinear ⇒ `eq_circumcenter_of_equidistant`; φ* = 0 ⇒ pole of a pair | `IsCompact.exists_isMinOn` (`Topology/Order/Compact.lean`), `isCompact_sphere`; repo `Developability/MaxCovariance.lean` (`circumcenter_equidistant`, `eq_circumcenter_of_equidistant`) |
| dev-min-width-hull-edge | Houle–Toussaint 1988 | partial (intrinsic max is the min) | **A-M** (~400–700) | statement the code needs: min over unit u of width is attained at u ⊥ (P_b − P_a) for some pair a<b (the code enumerates *all* pairs); proof: at a min, if max and min are each attained once, width is locally R·cos(t−t₀) > 0, no local min; no convex hull needed | `IsCompact.exists_isMinOn`; repo `Developability/Width.lean` (`width_exists_maxTo1`, `width_exists_minTo1`, `pairWidthSq_eq_sq_widthAlong`) |
| dev-lbfgs-two-loop | N&W Alg. 7.4 | N (efficiency; runtime gᵀd<0 check) | **A-S, already proved** | instantiate `RS.twoLoop_eq` and `RS.bfgsInv_posDef` with H₀ = γI, γ = sᵀy/yᵀy > 0; tie the iterative Rust loop (oldest-first deque) to the recursive `twoLoop` (newest-first list) | repo `RepulsiveSurfaces/Optimize.lean` (`twoLoop`, `bfgsInv`, `twoLoop_eq`, `bfgsInv_posDef`); `Matrix.PosDef` |
| dev-lewis-overton-termination | L&O 2013 | N (search capped; outcomes proved) | C as stated (semi-algebraic, Clarke) | A-M: for C¹ φ along the ray, existence of an Armijo + weak-Wolfe step (MVT) and finite termination of doubling/bisection (limit argument with continuity of φ′) | `exists_hasDerivAt_eq_slope`; Rademacher now exists (`LipschitzWith.ae_differentiableAt`, `Analysis/Calculus/Rademacher.lean`) for the Lipschitz version; repo `Developability/Optimizer.lean` (`lineSearchLoop_ok`) |
| dev-armijo-backtracking | N&W Alg. 3.1 | N (capped by `max_shrinks`) | A-S (<100) | smooth case; more useful: f = max of C¹ pieces, Armijo eventually holds iff gᵢᵀd < 0 for every *active* piece (shows the chosen-piece gradient is not enough) | `HasDerivAt.tendsto_slope…`, `tendsto_pow_atTop_nhds_zero_of_lt_one` |
| dev-flip-manifold | Edelsbrunner 2001 | N (mesh is re-validated) | A/B-M on the repo's `TriMesh` model | flip keeps |V|, |E|, |F| (so χ) and edge degree ≤ 2 when k≠l and kl is not an edge | repo `GaussBonnet.lean` (`faceEdges`, `edges`, `edgeDegree`, `eulerChar`, `Manifold`) ; **Status (2026-11-05): trimmed** — the topological part proved on gluing data (`Delaunay/Complex.lean`, `Delaunay/Flip.lean`, `it-complex-flip`: V with classes, E, F, χ, boundary, one fan per vertex kept; any coincidences of i, j, k, l) and transferred to oriented fan-connected `TriMesh` (`Delaunay/FromTriMesh.lean`); still cited: the flipped mesh is simplicial and `Manifold` (where k≠l, kl∉E enter) |
| dev-link-condition | Dey+ 1999 | N (re-validated) | combinatorial part B-L; homeomorphism C | link condition ⇒ contracted complex has cycle/path vertex links | Mathlib `AlgebraicTopology/SimplicialComplex/Basic.lean` (`AbstractSimplicialComplex`, `map` — contraction is a vertex `map`); no links / combinatorial manifolds |
| dev-cut-edge-weight | Dev. Sec. 4.4 | N (design choice) | — (nonnegativity is one Cauchy–Schwarz line) | not worth a theorem | `abs_real_inner_le_norm` |

#### Shared infrastructure worth building once

**Integral over a triangle** (B-M, ~600–1000 lines): affine change of variables to the standard simplex
(`integral_image_eq_integral_abs_det_fderiv_smul`, or `addHaar_image_linearMap` for the area), Fubini over
{x,y ≥ 0, x+y ≤ 1} (`setIntegral_prod`), the Beta integral for monomials. It would cover the area-coordinate rule
(app-vector-fem-connection-laplacian), the identification "these entries are ∫φᵢφⱼ" left out of
`app-fem-galerkin-mass` and `cotan-fem-stiffness`, Chen–Holst Lemma 4.1 (it-odt-energy), and the L² norms in
the AMR estimator. **Status (2026-11-05): built** (`DdgProofs/TriangleIntegral/`; per-item status in the
"Shared infrastructure" list near the top). Route taken: Fubini via `integral_prod` on the indicator, the Jacobian
theorem with `det = sda`, the Beta integral by induction on the exponent (no Gamma functions needed).

#### Notes (only where not obvious)

- **dev-simple-eigenvalue-derivative.** The Rust tag sits on `covariance` (`energy/covariance.rs:21`), which uses
  only λ_min. The repo's `IsEigenbasis` for `Sym3` gives the expansion u = Σ⟨u,eᵢ⟩eᵢ. For unit u = c·x₀ + s·z
  (z ⊥ x₀, c²+s² = 1) we get uᵀ(A+H)u ≥ λ + x₀ᵀHx₀ + g s² − 2‖H‖|s| − 2‖H‖s², so
  λ_min(A+H) − λ − x₀ᵀHx₀ ≥ −‖H‖²/(g − 2‖H‖). With `eigenvalue_supergradient` as the upper bound this gives
  `HasDerivAt` along any C¹ curve A(t). This is the only entry in this share where a gradient formula the code
  uses depends on a cited result.
- **dev-min-width-hull-edge / dev-max-covariance-voronoi.** In both, the code enumerates a superset of the
  paper's candidates (all pairs; all signed triples plus poles). So the statement to prove is "some global
  minimizer lies in the candidate set". Proving that needs only compactness, finite max/min and
  one-dimensional calculus, not hulls or Voronoi diagrams. Width is the easier of the two (circle, 1-D);
  do it first and reuse the active-set argument for the sphere.
- **app-steiner-fermat-point.** The minimality of the point is the easy part: Σ‖x−pᵢ‖ ≥ Σ⟨x−pᵢ, uᵢ⟩ = Σ‖F−pᵢ‖
  when Σuᵢ = 0. The work is the algebra showing that the code's barycentrics give 120° angles. An easier route
  is to state it as dⱼ² + dₖ² + dⱼdₖ = aᵢ² (law of cosines at 120°) with the d's from the proved
  `cornerDistances`.
- **app-vector-singularity-index.** The Rust doc (`vector.rs:478`) uses ω ∈ (−π, π]; CITED says (−π, π). With
  the half-open interval, an edge with ω_ij = π has ω_ji = π as well, so the twin terms do not cancel. The
  Poincaré–Hopf sum nχ then needs "no edge at exactly π". A formal proof would have to state this hypothesis.
- **it-straightest-geodesic.** Most of the multi-step statement already exists in pieces (one step in
  `Trace.lean`, rigid motions in `Diamond.lean`, the strip crossing in `IntegerCoords/Strip.lean`). The new part
  is an inductive "strip" record (a list of faces given by lengths, plus the shared edges) and the lemma that the
  ray stays one line under the composed unfoldings. The halfedge bookkeeping stays unverified.
- **it-rippa-theorem.** Proof route: E₁(f+a) − E₁(f) − E₁(a) = 2∇a·∮f n, and the right side is the same for
  both triangulations, so the difference Q = E₁ − E₂ is a quadratic form vanishing on affine data with
  B(f,a) = 0. Hence Q has rank ≤ 1, and its sign is computed at f = e_k.

#### Top candidates (value to library correctness × feasibility)

1. **dev-lbfgs-two-loop**: S, already proved as `rs-lbfgs-two-loop`; needs only an instance with H₀ = γI and the tag.
   (The same holds for `rc-lbfgs-two-loop`, CITED.md ~l.1841, outside this share.)
2. **app-amr-edge-flux**: S; a library formula (AMR jump residuals); it is `divergence_corner` × 2.
3. **dev-simple-eigenvalue-derivative**: M; the gradient of E_λ the code uses; infrastructure mostly in `Covariance.lean`.
4. **dev-min-width-hull-edge**: M; turns "upper bound" into "equals the minimum" for the intrinsic max energy.
5. **dev-max-covariance-voronoi**: M-L; the same for λ_max; reuses `eq_circumcenter_of_equidistant`.
6. **app-steiner-fermat-point**: S for minimality, M for the 120° / S² algebra; place and length of the Steiner points.
7. **it-straightest-geodesic** (multi-step strip): M; the tracing query as a whole (currently only one step).
8. **Triangle-integral infrastructure**: M; covers four cited "∫" gaps (FEM mass, connection Laplacian, ODT Lemma 4.1, cotan stiffness).
9. **app-vector-singularity-index**: per-face integrality S; Poincaré–Hopf M (makes the π-tie hypothesis explicit).
10. **it-rippa-theorem** (one diamond) and **it-delaunay-flip-termination** (one-flip potential decrease): M each, low dependency.

Cheap but low value: app-vector-maximum-principle, dev-armijo-backtracking, app-amr-bulk-marking (greedy),
dev-valence-three-flat (after A.1). Not worth it: it-odt-target-length-splits, dev-cut-edge-weight,
app-vector-transport (definitions or heuristics). Out of reach: it-delaunay-existence-uniqueness,
Chew/Ruppert termination, the AMR estimator bounds and contraction, Chen–Holst Thm 2.3, the homeomorphism part
of dev-link-condition.

#### Stale or inaccurate CITED.md reasons

1. **dev-lbfgs-two-loop** (l.1416): "linear algebra of rank-two updates" is given as the reason, but it is
   already machine-checked generically (`DdgProofs/RepulsiveSurfaces/Optimize.lean`: `twoLoop_eq`,
   `bfgsInv_posDef`, any SPD H₀). The entry should cross-reference it or be proved by instantiation.
   (Same for `rc-lbfgs-two-loop`.)
2. **app-amr-edge-flux** (l.1228): "plane geometry of one triangle; cotan-fem-stiffness is cited too". The cited part
   of `cotan-fem-stiffness` (Lebesgue integral, convergence) is not needed. The identity is `divergence_corner`
   (HeatFem.lean) with X = ∇u, and `cotan_formula` is proved.
3. **app-steiner-fermat-point** (l.1129): "nothing in lean/proofs is set up for the trigonometry of the triangle's
   angles" is out of date. `Triangle.lean` has `cornerAngle` (arccos), `cos_cornerAngle`, `sin_cornerAngle`,
   `cornerAngle_add_cornerAngle_add_cornerAngle`, and `angle_eq_cornerAngle` linking to Mathlib's
   `EuclideanGeometry.angle`. `Diamond.lean` has `diamondX_eq_cos` and `diamondH_eq_sin`.
4. **dev-min-width-hull-edge** (l.1399): "needs convex hulls of planar point sets" is overstated. The code's
   candidate set is all pairs, so the required statement has no hull in it.
5. **dev-max-covariance-voronoi** (l.1381): "would need the cell decomposition of the sphere" is overstated. An
   active-set argument at a minimizer, with compactness (`IsCompact.exists_isMinOn`), suffices.
6. **dev-simple-eigenvalue-derivative** (l.1358):
   - "Mathlib has … the min-max (Rayleigh) characterization": Mathlib has only the extremal case
     (`hasEigenvalue_iInf/iSup_of_finiteDimensional`, `Rayleigh.lean`). There is no Courant–Fischer (grep
     "courant|fischer|minmax" finds nothing).
   - The conclusion "no perturbation theory, so not formalizable" does not hold for the smallest eigenvalue,
     which is the only case the code uses.
7. **app-vector-singularity-index** (l.1202):
   - The interval (−π, π) disagrees with the code's (−π, π] (`vector.rs` doc).
   - "needs … Gauss–Bonnet": Gauss–Bonnet is already proved (`GaussBonnet.lean`, `gauss_bonnet_closed`). What
     is missing is an orientation hypothesis on `TriMesh`.
8. **dev-flip-manifold** (l.1463): "combinatorics of the mesh data structure, which lean/proofs does not model" is
   partly stale. `GaussBonnet.lean` models faces, `Sym2` edges, edge degree, boundary and χ, which is enough for
   "a flip keeps χ and edge degree ≤ 2". It does not model vertex links.
9. **dev-link-condition** (l.1477): the reason is still accurate, but Mathlib now has `AbstractSimplicialComplex`
   with `map` (`AlgebraicTopology/SimplicialComplex/Basic.lean`); it has no links and no notion of
   combinatorial manifold.
10. **dev-lewis-overton-termination** (l.1434): the reason is accurate for the semi-algebraic theorem. Note that
    Rademacher is in Mathlib (`LipschitzWith.ae_differentiableAt`), and that the C¹ special case needs no
    nonsmooth analysis.
11. **app-vector-fem-connection-laplacian** (l.1194): "Mathlib has the integrals" is loosely true. There is no
    integral over a triangle; it needs change of variables plus Fubini (infrastructure above). The code's actual
    formula (the area-coordinate rule) is formalizable without the gauge model.

Accurate as written: it-delaunay-existence-uniqueness, it-delaunay-flip-termination (surface), it-chew, it-ruppert
(no Delaunay in Mathlib at all), it-rippa ("not asked for"), it-odt-*, app-amr-residual-estimator (no domain
Sobolev spaces or Clément interpolation; Mathlib Sobolev is Fourier/Bessel only), app-amr-bulk-marking,
dev-hinge-* (embeddedness genuinely absent, though algebraic stand-ins suffice), dev-armijo-backtracking,
dev-cut-edge-weight, app-vector-transport, app-vector-maximum-principle.

### CITED.md survey, part 3 (lines 1495–2630)

Sections: Repulsive Curves (RC), Integer Coordinates (IC), Intrinsic Simplification (IEM),
Repulsive Surfaces (RS), Points as Tori (PAT), Out of scope. Pinned Mathlib = `lean/proofs/.lake/packages/mathlib`
(v4.34.1). Every Mathlib name below was found by grep, and the starred (★) ones were also `#check`ed in a
scratch file (`lake env lean`, nice 19; all resolved).

Depends: **Y** = a library output or safety claim depends on it; **P** = partly (an output
that is reported or tested, or one branch); **N** = motivation, efficiency, a baseline, or test data only.
Class: **A** = provable now; **B** = needs bounded missing infrastructure; **C** = out of reach. Effort:
S < 300 lines, M 300–1500, L > 1500.

#### Table

| id | paper | dep | class/effort | partial option | key Mathlib / DdgProofs |
|---|---|---|---|---|---|
| rc-exponent-range | RC A.2 / Blatt | N | C | none needed (`Exponents` bookkeeping is already the only use) | no fractional Sobolev spaces |
| rc-differential-order | RC A.3 | N | C | none | — |
| rc-fractional-operator | RC A.1/A.4 | N | C | none | — |
| rc-quadrature-consistency | RC 5.1 | N | C | B-L: omitted adjacent-pair terms → 0 for one fixed smooth family (circle); low value | — |
| rc-chord-convergence | Ortega–Rheinboldt / Kelley | P (failures rejected) | B-L | A-M: one-step residual bound `‖Φ(γ̃+x)‖ ≤ ½L‖x‖² + ‖C(γ̃)−C(γ)‖‖x‖`; contraction then gives local linear convergence | ★`Convex.norm_image_sub_le_of_norm_fderiv_le`, ★`ContractingWith.efixedPoint`, ★`HasStrictFDerivAt.implicitFunction`; DdgProofs `rc-chord-projection` |
| rc-distance-gradient | Bauschke–Combettes 12 | **Y** (∇ of surface potential; the last missing link of `rc-penalty-surface`) | **A-S** (~200–300) | full statement for any closed convex `K` in a real inner-product space | ★`norm_eq_iInf_iff_real_inner_le_zero`, ★`exists_norm_eq_iInf_of_complete_convex`, ★`hasStrictFDerivAt_norm_sq`, ★`HasFDerivAt.sqrt`; `RC.penalty_surface` |
| rc-closest-points | Ericson 5.1.5/5.1.9 | **Y** (CCD safety: `primDist_advance` is proved for the true infimum `primDist`; the code uses these formulas) | **A-M** (600–1000) | each `ccd.rs` distance = `primDist` (exact reals); segment–segment exact except when `0 < denom ≤ 1e-14·a·e`, where the overestimate is bounded (see notes); obstacle Ericson routine for nondegenerate triangles | ★`norm_eq_iInf_iff_real_inner_le_zero` (or direct convex-quadratic algebra); `RC.primDist`, `RC.Weights` (Ccd.lean) |
| rc-near-cos-series | DLMF 8.17/15.2 | **Y** (value of the near-pair inner integral `F(y)`) | **A-M** (500–900) | series identity plus an explicit bound on the 60-term remainder; plugs into the existing `near_inner_integral` | ★`Real.one_add_rpow_hasFPowerSeriesOnBall_zero`, ★`Real.one_div_one_sub_rpow_hasFPowerSeriesOnBall_zero`, ★`hasDerivAt_tsum`, ★`intervalIntegral.hasSum_intervalIntegral_of_summable_norm`, ★`Ring.choose` |
| rc-near-closest | Ericson; Danskin | P (gradient of the blended energy; tested with finite differences) | A-M (fixed active set) / B-L (Danskin at switch points) | interior case: explicit rational `s*,t*`, the envelope `∇d = n̂⊗(…)` from stationarity; clamped cases with the active set locally constant | ★`IsLocalMin.hasFDerivAt_eq_zero`; no Danskin/envelope theorem in Mathlib |
| rc-near-quadrature-error | Trefethen 19.3 | N (accuracy is measured) | C | trivial: `sinh(±iπ/2) = ±i` (pole location) | no Gauss–Legendre and no Bernstein-ellipse bounds; only Chebyshev–Gauss (`Polynomial.Chebyshev.integral_eq_sumZeroes`) |
| rc-armijo | N&W 3.1 | N (halvings capped) | **A-S** (~100) | also closes `dev-armijo-backtracking` (other share) | ★`hasDerivAt_iff_tendsto_slope`, ★`exists_pow_lt_of_lt_one` |
| rc-accd | Li et al. 2021 §5 | N (safety proved; cap of 100) | A-S | termination with the code's slack `1e-9·d₀`: at most `(b−a)l/slack` steps; counterexample with slack 0 (see stale list) | `RC.advance_safe` |
| rc-isotopy | RC 5.4; Burde–Zieschang | P (knot-type claim) | C | A-S: "embedded at every time" if nonadjacent pairs are CCD-checked **and** adjacent edges do not fold (a hypothesis the code does not check, as CITED says) | no knots or isotopy in Mathlib (only `Quandle`) |
| rc-bh-admissibility | RC B.1.2 | N (accuracy; θ→0 measured) | B-L | A-M: first-order bound `|Σ_B − lumped| ≤ sup‖∇k‖·Σ ℓ_j(‖x_j−x̄‖+‖T_j−T̄‖)` on the convex hull | ★`Convex.norm_image_sub_le_of_norm_fderiv_le`; `RC` `rc-bvh-lumped` |
| rc-bct-admissibility | RC B.2; Hackbusch | N (preconditioner quality) | B-L | the same mean-value bound per block | same |
| rc-multigrid | RC B.3; Notay | N (speed) | C | none worth doing | — |
| rc-lbfgs-two-loop | Liu–Nocedal; N&W 7.4 | N (efficiency) | **A-S** (~100–200; most of it already exists) | reuse `RS.twoLoop_eq`, `RS.bfgsInv_posDef`; add a PSD or subspace variant for `H₀ = γM⁻¹` (mean-free); also closes `dev-lbfgs-two-loop` | DdgProofs/RepulsiveSurfaces/Optimize.lean; `RC.massRows_posDefOn` |
| rc-ncg-polak-ribiere | PR 1969; N&W 5.2 | N | A-S (descent via restart, trivial) / B-L (global convergence) | only the restart rule is worth doing; the convergence theorems need strong Wolfe conditions, and the code uses Armijo only | — |
| ic-normal-isotopy | Kneser/Haken | **Y** (completeness of the encoding) | C | B-L, combinatorial: existence and uniqueness of a whole-complex model for valid `n` on the existing `Cx` model (Tier C) | no Jordan–Schoenflies, Δ-complexes or isotopy (grep confirms) |
| ic-refinement-termination | IC Thm C.1 | N (insertion limit) | C | none worth doing (planar Shewchuk packing alone is L+) | no Delaunay/Voronoi in Mathlib |
| ic-dijkstra-ball | Xia 2013 | N | C | none (weaker stretch bounds are also long) | — |
| ic-mollify | Sharp–Crane 2020 | N | n/a | the formula is already proved (`intrinsic-mollification`); what remains is an empirical claim | — |
| ic-transfer-mass-integral | Ciarlet 4.1 | N (finite identity proved; Rust checks it) | B-M | `∫_T λ_aλ_b = A(1+δ_ab)/12` via Fubini and an affine change of variables; additivity over a.e.-disjoint convex pieces (tiling as a hypothesis) | ★`MeasureTheory.integral_image_eq_integral_abs_det_fderiv_smul`, ★`Convex.addHaar_frontier`, `MeasureTheory.integral_prod` |
| ic-vertex-removal-delta | IC App. D (open) | P (refusals) | open | instance checks: whether the stuck stars of SPEC §4.6 (vertex with a self-edge) admit a bounded-length flip sequence, decided on small `Cx` models; could yield a counterexample | `Cx.flip_valid`, `FlatStar.exists_flippable` |
| iem-wasserstein-memory | IEM 5.2 | N | n/a (no bound is claimed) | A-S/M: memory cost ≤ cost of the merge-history coupling (Jensen) and W₁ ≤ the same; counterexamples showing that neither of memory cost and W₁ bounds the other | `IsCoupling`, `wasserstein_eq`, `centerOfMass_combine`; ★`ConvexOn.map_sum_le` |
| iem-cetm-existence | CETM 3.1 | N (only the one-vertex case is used, and it is proved) | B-L | M: monotonicity `⟨Θ(u)−Θ(v),u−v⟩ ≤ 0` from a cotan-Laplacian Jacobian on the regular region; low value | no Lobachevsky/Clausen function (grep: only von Staudt–Clausen); multivariate convexity only by restricting to lines |
| iem-delaunay-quality | IEM 6.3 = IC C.1 | N | C | as `ic-refinement-termination` | — |
| rs-exponent-range | RS 2.1 | N | C | none needed | — |
| rs-differential-order | RS 2.3 | N | C | none | — |
| rs-fractional-operator | RS 2.4 | N | C | none | — |
| rs-preconditioner-approx | RS 5.2 | N (GMRES iteration count) | C | none (the exact matrix identity is proved) | — |
| rs-preconditioner-scale | Di Nezza et al. Prop 3.3 | N (iteration count) | B-L | B-M: the radial–angular reduction by polar coordinates; the Mellin transform of `1−cos` is missing (`Complex.integral_cpow_mul_exp_neg_mul_Ioi` needs `Re a > 0`, and `a = i` is the boundary case) | ★`Real.Gamma_mul_Gamma_one_sub`, ★`integral_comp_polarCoord_symm`, `integral_sin_pow` |
| rs-aqp | Kovalsky et al. | N (comparison only) | n/a (heuristic) | none | — |
| rs-bcqn | Zhu et al. | N (comparison only) | n/a | none (L-BFGS part is proved) | — |
| rs-energy-consistency | RS 8.1 | N | C (empirical) | none | — |
| rs-chen-holst | Chen–Holst 2011 | N | n/a (heuristic) | A-M: the smoothing target is the stationary point of the local ODT energy; low value | — |
| pat-network | PAT 4.3 | N | n/a (a definition) | — | — |
| pat-laplace-method | PAT 2.2 | N (motivation) | C (the continuous version) | **A-S**, discrete Eq. 25 as evaluated: `|φ − g_{i*}| ≤ 2·max|g|·Σ_{i≠i*}(a_i/a_{i*})e^{−λ(ρ_i−ρ_{i*})}` and `−λ⁻¹log Σ a_i e^{−λρ_i} → min ρ` | `PAT.blend`, `PAT.blend_convex`; `Real.exp_le_exp`, `Real.log_le_log` |
| pat-signed-hopf-cole | PAT 2.1, App. A–B | N | C | none | — |
| pat-winding-dipole | Barill 2018 | N (baseline, sign of the ground truth) | C | none beyond the dipole derivatives already proved | — |
| pat-point-area | Barill 2018 §5 | N | n/a | — | — |
| pat-truth-closest-point | Ericson 5.1.5 | P (ground truth of the measurements) | **A-M**, shared with rc-closest-points | the same Ericson routine (the PAT copy has degenerate guards) | as rc-closest-points |
| pat-solid-angle | van Oosterom–Strackee | P (truth sign; FWN triangle elements) | B-M/L | M: define the signed excess with dihedral angles and prove the van Oosterom `atan2` form; also closes `spherical-triangle-area` (other share). Girard's theorem as a sphere measure is L; winding-number integrality C | ★`InnerProductGeometry.angle`, ★`Complex.arg`; no spherical trigonometry (grep: no Girard/solidAngle) |
| pat-triangle-sampling | Osada et al. | N (training distribution) | B-M | A-S core: the Jacobian of `(u,v) ↦ √u(1−v)a+√u v b+(1−√u)c` has constant absolute value = Area | ★`MeasureTheory.integral_image_eq_integral_abs_det_fderiv_smul` |
| pat-query-sampling | reference code | N | n/a (a definition) | — | — |
| pat-regularized-wn | Chen et al. 2024 | N (baseline) | B-M (no `erf` in Mathlib) | A-S: with erf defined as `∫₀ᵗ`, `S′(t) = (4/√π)t²e^{−t²}` (FTC), hence `S = O(t³)` and the kernel is smooth at ρ = 0 | `integral_gaussian`, FTC |
| pat-lse | Madan–Levin | N (baseline) | **A-S** (<100) | the full bound, over the neighbor set (code sums over `nbrs`, so `log n` = log of the neighbor count) | `Real.log_le_log`, `Finset.single_le_sum` |
| pat-adam | Kingma–Ba | N | n/a | — | — |
| pat-offset | Ricci; Quilez | N (apps claim no exactness) | B-M (no signed distance in Mathlib) | A-S: unsigned `infDist x (A∪B) = min`, offset level set ⊆ frontier of `cthickening` | ★`Metric.infDist`, ★`Metric.frontier_cthickening_subset`, `Metric.infDist_le_infDist_of_subset` |
| pat-quilez-2d | Quilez 2019 | N (test ground truth) | A, S each for segment/box/pie/arc; M–L for all 16 (ellipse: cubic solve) | follow the `PAT.circle_sdf` pattern ("≤ distance to every point, attained") | Pat2d.lean |
| (xref) pat-fwn-expansion error | Barill §4 | N | B-L | M: order-0 bound from the mean-value inequality on `1/(4π|y|)` | `Convex.norm_image_sub_le_of_norm_fderiv_le`, `taylor_mean_remainder_lagrange` |
| (xref) pat-torus-sdf spindle | Quilez 2025 | N | A-M | Eq. 23 = exact distance to the tube of the generating circle for `R < |r|` | Torus.lean |
| pat-render-ggx / pat-render-fresnel | demo shader | N | n/a | — | — |

#### Notes (only where non-obvious)

- **rc-closest-points.** `ccd.rs::point_segment_distance` clamps the projection, which gives the exact
  distance. `point_triangle_distance` returns the minimum of the three edge distances and, if the projection lies inside,
  the plane distance. That is exact for every triangle, collinear ones included (their hull is the union of their
  edges). `segment_segment_distance` returns the minimum of the 4 endpoint distances and the interior pair, but
  only when `denom > 1e-14·a·e`. Every candidate is the distance of a real pair, so the result
  ≥ true distance, which is the unsafe direction for CCD. Moving from the interior minimizer
  along the valley `t(s)` to the boundary of `[0,1]²` raises the squared distance by `(denom/e)(Δs)² ≤ denom/e`.
  So the expected theorem is `cand² ≤ d² + 1e-14·|q₁−p₁|²` in the skipped window. This is the precise form
  of the code comment "second-order amounts". Together with `primDist_advance`, the theorem gives an
  end-to-end exact-real CCD safety statement, with a quantified tolerance. `obstacle.rs` and PAT
  `truth.rs` each carry their own copy of Ericson's Voronoi routine. For a nondegenerate triangle the final branch has
  `va+vb+vc = |ab×ac|²` (Lagrange identity), and the edge denominators are `|ab|²`, `|ac|²`, `|bc|²`. The RC copy has
  no degenerate guards, but `PreparedObstacle::new` rejects degenerate triangles, so that is fine.
  **Status (2026-11-05): done** (`ClosestPoints.lean`): exact distances = `primDist`, end-to-end CCD safety
  (`segSeg_advance_safe`, `ptTri_advance_safe`), the Voronoi routine of `truth.rs` for every triangle
  (`closestPtTri_isClosest`, via the variational inequality and region coverage), the `obstacle.rs` copy equal to it
  for non-degenerate triangles (`closestPtTriRC_eq`).
- **rc-distance-gradient.** The proof is short. Let `f = ½‖x−π x‖²`. Upper bound: `f(x+h) ≤ ½‖x+h−πx‖²`.
  Lower bound: `f(x) ≤ ½‖x−π(x+h)‖²` together with nonexpansiveness of `π`, which follows from the variational inequality. Both sides
  are `f(x)+⟨x−πx,h⟩ ± O(‖h‖²)`. Then `δ = √(2f)` follows by `HasFDerivAt.sqrt`. Mathlib has **no convex-projection map**,
  only existence and the characterization, so the proof takes `π` as any function satisfying the VI.
- **rc-near-cos-series.** The derivative of `Σ C(a,k)(−1)^k s^{2k+1}/(2k+1)` is the binomial series at `−s²`, which
  Mathlib now has for real exponents with radius ≥ 1. The second series is
  `one_div_one_sub_rpow_hasFPowerSeriesOnBall_zero` with `a = ½` exactly. Composing with `sin` and applying the existing
  `near_inner_integral` (any `S` with `S′ = cos^m`) gives the identity. The truncation bound is elementary
  (`|C(a,k+1)/C(a,k)| = |a−k|/(k+1)`, `x ≤ ½`).
- **rc-lbfgs-two-loop.** RC uses `H₀ = γ_k M_k⁻¹` (mean-free solve), which is PSD and PD only on `ker M_c`. The RS proof
  needs PD on all of ℝⁿ, so a 30–80-line variant is needed: PSD in gives PSD out, and PD on an invariant subspace
  containing all `s_i, y_i`. The same file closes `dev-lbfgs-two-loop` (`H₀ = γI`).
- **rc-armijo.** Only `HasDerivAt φ φ′ 0` with `φ′ < 0` is needed (C¹ is not). The value is low: the projected and
  CCD-capped step is not `γ+τd`, and the halvings are capped.
- **iem-wasserstein-memory.** Small counterexamples (1D, two coarse vertices with masses assigned crosswise, or a symmetric merge) show
  that neither `Σ m_j‖t_j‖ ≤ W₁` nor `≥` holds. The only honest inequalities are the two "≤ history-coupling cost"
  bounds. They are cheap with the existing finite-coupling model and of low value.
- **pat-laplace-method / pat-lse.** One lemma about `−λ⁻¹ log Σ a_i e^{−λρ_i}` and the self-normalized average covers both.
  It is the only exact statement here about Eq. 25 as the library evaluates it in the large-λ regime.
- **Out of scope section.** Nothing listed there is a mathematical claim. The only items in reach are
  re-modelled algorithms with no mechanical tie to Rust, e.g. "the BFS satisfies `CoarsenRule`"
  (RC header) or the iteration and stop-rule bounds. They are A-S/M each and of low value. The render entries make no claim.

#### Top candidates (value to library correctness × feasibility)

1. **rc-closest-points (+ pat-truth-closest-point)**: A-M. Joins the proved CCD advancement
   (`primDist`) to the distances the code computes. It is the only safety claim in this share with a feasible proof, and the same
   lemmas give the PAT ground truth.
2. **rc-distance-gradient**: A-S. Completes the surface-potential gradient chain (`rc-penalty-surface`
   is proved "given ∇δ").
3. **rc-near-cos-series**: A-M. The value of every near-pair inner integral. The Mathlib support is now
   there (the CITED reason is stale).
4. **rc-lbfgs-two-loop (+ dev-lbfgs-two-loop)**: A-S. Most of it is already proved in `RS/Optimize.lean`.
5. **rc-near-closest (fixed-active-set case)**: A-M. The gradient of the blended near-pair energy away from the
   active-set switches.
6. **rc-accd termination + slack counterexample**: A-S. Corrects a misstatement (below).
7. **rc-armijo (+ dev-armijo-backtracking)**: A-S. Trivial and gives a finite-halving bound.
8. **pat-lse / discrete pat-laplace-method**: A-S. Cheap and exact about Eq. 25.
9. **pat-solid-angle (+ spherical-triangle-area)**: B-M. Covers two ids, but the trigonometry with `atan2` branches is fiddly.

#### Stale or incorrect CITED.md reasons

- **rc-near-cos-series**: "binomial series with real exponents are in Mathlib only in pieces" is stale.
  `Real.one_add_rpow_hasFPowerSeriesOnBall_zero` (radius 1) and `Real.one_div_one_sub_rpow_hasFPowerSeriesOnBall_zero`
  (Analysis/Analytic/Binomial.lean) give both series. Termwise work goes through `hasDerivAt_tsum` and
  `intervalIntegral.hasSum_intervalIntegral_of_summable_norm`.
- **rc-lbfgs-two-loop**: "As for `dev-lbfgs-two-loop` (rank-two updates)" is stale. The rank-two algebra and
  positive definiteness are machine-checked as `rs-lbfgs-two-loop` (`RS.twoLoop_eq`, `RS.bfgsInv_posDef`). Only a
  subspace/PSD variant is missing. The same staleness applies to `dev-lbfgs-two-loop` (other share).
- **rc-accd**: "in finitely many advancements for a positive gap" is false for the crate's rule
  `τ ← τ + (d − g)/l` (the rule of `RC.advance_safe`) in exact arithmetic. With `d(τ) = g + l′(τ*−τ)`, `l′ < l`, the
  iterates converge geometrically to `τ*` and never stop. Termination comes from the slack `1e-9·d₀` in
  `ccd.rs::advance` (at most `(b−a)l/slack` steps, about 10⁹, hence the cap of 100), not from the gap.
  Li et al.'s ACCD steps by a fraction of `d` and stops below the gap; that variant does terminate.
- **rc-distance-gradient**: "Mathlib has the orthogonal projection onto … convex sets" is imprecise. Mathlib has
  orthogonal projection only onto subspaces (`starProjection`). For convex sets it has only the existence
  (`exists_norm_eq_iInf_of_complete_convex`) and the characterization (`norm_eq_iInf_iff_real_inner_le_zero`), with
  no projection map. The conclusion (not formalized) stands, but it is A-S.
- **rc-armijo / dev-armijo-backtracking**: the stated hypothesis "continuously differentiable near γ" is stronger than needed:
  differentiability at γ is enough.
- **rc-ncg-polak-ribiere**: "converges globally under the usual line-search conditions" is not wrong, but those conditions
  (strong Wolfe and sufficient descent, Gilbert–Nocedal 1992 for PR+) are not met by the crate's Armijo-only search.
  The cited references (PR 1969, N&W 5.2) are not the PR+ convergence source.
- All other reasons in this share check out against the pinned Mathlib. Searches found no fractional Sobolev
  spaces, Wasserstein distance, Lobachevsky/Clausen function, erf, Gauss–Legendre/Bernstein-ellipse bounds, Danskin
  theorem, spherical trigonometry/Girard, Jordan–Schoenflies, isotopy/knots, Delaunay/Voronoi, or signed distance.
