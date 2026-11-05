import DDG.Repulsive.Curves.BarnesHutAdjoint
import DDG.Repulsive.Curves.Constraints
import DDG.Repulsive.Curves.DualLengths
import DDG.Repulsive.Curves.Energy
import DDG.Repulsive.Curves.Gram
import DDG.Repulsive.Curves.Hierarchy
import DDG.Repulsive.Curves.Kernel
import DDG.Repulsive.Curves.KernelGradient
import DDG.Repulsive.Curves.Multigrid
import DDG.Repulsive.Curves.Penalties
import DDG.Repulsive.Curves.Saddle
import DDG.Repulsive.Surfaces.BarnesHut
import DDG.Repulsive.Surfaces.Basic
import DDG.Repulsive.Surfaces.Constraints
import DDG.Repulsive.Surfaces.Energy
import DDG.Repulsive.Surfaces.Geometry
import DDG.Repulsive.Surfaces.Kernel
import DDG.Repulsive.Surfaces.Metric
import DDG.Repulsive.Surfaces.Penalties
import DDG.Repulsive.Surfaces.Precond
import DDG.Repulsive.Surfaces.Saddle

/-!
# DDG.Repulsive: tangent-point energies of curves and surfaces

**Scope.** Machine-checked derivations of the formulas of Repulsive Curves and Repulsive Surfaces:
the discrete tangent-point energies and their exact differentials along any motion, the fractional
Sobolev inner products, saddle-point systems, constraints and penalties, Barnes–Hut and the
multigrid hierarchy. Curves are in namespace `DDG.RC`, surfaces in `DDG.RS`; vectors and the
curve-motion calculus are `DDG.Geometry.Vec3`. The smooth theory (Blatt's lemma, the order of the
operators), the error bounds of Barnes–Hut and of the block cluster tree and all convergence
statements are cited, not proved (`lean/proofs/CITED.md`).

**Main results.**
* `Curves/Kernel`, `KernelGradient`, `Energy`: the kernel `k^α_β` (paper Eq. 3), the discrete
  energy (Eq. 17–18) and its differential (`hasDerivAt_pairEnergy`, `energyGrad`), invariances.
* `Curves/Gram`: assembly and positive semidefiniteness of the Sobolev–Slobodeckij Gram matrices
  (`gramB_posSemidef`). `Saddle`: saddle systems, the projector and the scaling of the direction.
* `Curves/Constraints`, `Penalties`, `DualLengths`; `Hierarchy`, `BarnesHutAdjoint`,
  `Multigrid` (the black/white coarsening rule `no_adjacent_whites`, prolongation).
* `Surfaces/*`: the same for surfaces: kernel, energy and its gradient
  (`energy_gradient_identities`), the operators of the inner product (`Metric`), the fractional
  preconditioner order (`precond_order`), saddle systems and the Schur complement, constraints,
  penalties, the Barnes–Hut energy (`bh_energy_split`).

**Sources.** C. Yu, H. Schumacher, K. Crane, *Repulsive Curves*, ACM TOG 40(2), 2021;
C. Yu, K. Brakensiek, H. Schumacher, K. Crane, *Repulsive Surfaces*, ACM TOG 40(6), Article 268,
2021. The Rust crates were developed with the authors' MIT-licensed implementations
(`repulsive-curves`, `repulsive-surfaces`) as reference (see the crate READMEs); the Lean statements
are written from the papers' text and the contracts in `crates/ddg-repulsive-curves/SPEC.md` and
`crates/ddg-repulsive-surfaces/SPEC.md`.

**Depends on:** `DDG.Geometry` (vectors, closest points and conservative advancement live there).
-/
