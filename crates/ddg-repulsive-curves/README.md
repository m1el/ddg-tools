# ddg-repulsive-curves

A Rust implementation of

> Chris Yu, Henrik Schumacher, Keenan Crane. *Repulsive Curves.* ACM Transactions on Graphics
> 40(2), Article 10, 2021. <https://doi.org/10.1145/3439429>

(with the paper's supplemental material): the tangent-point energy of curve networks in R³,
its descent in the fractional Sobolev–Slobodeckij metric of the paper, hard constraints and soft
penalties for design, collision-guarded time steps, and the Barnes–Hut, hierarchical-matrix and
multigrid acceleration of the paper's Appendix B.

**Status: contract only.** The specification is [`SPEC.md`](SPEC.md): every discrete formula
with its equation number and `rc-` formula id, the places where the paper or the reference code
is wrong or silent (decisions D1–D30), the API, the defaults, the verification and measurement
plans, and the work packages. The library is not implemented yet.

## Credits

* The method, the discretization, the metric, the constraints and penalties, the acceleration
  scheme, the knot datasets (Knot128, Trefoil100) and the Freedman–He–Wang meshes are by
  **Chris Yu, Henrik Schumacher and Keenan Crane**.
* The authors' reference implementation (`repulsive-curves`, MIT, Copyright (c) 2019 Christopher
  Yu) is read and may be ported; its licence notice is in
  [`LICENSES/repulsive-curves-MIT.txt`](LICENSES/repulsive-curves-MIT.txt). Code ported from one
  of its dependencies will carry that dependency's notice in `LICENSES/` as well.
* The tangent-point energy and its energy space: S. Blatt, *The energy spaces of the tangent
  point energies*, J. Topol. Anal. 5 (2013); S. Blatt, Ph. Reiter, *Regularity theory for
  tangent-point energies: the non-degenerate sub-critical case*, Adv. Calc. Var. 8 (2015).
* Conservative continuous collision detection: M. Li, D. M. Kaufman, C. Jiang, *Codimensional
  Incremental Potential Contact*, ACM TOG 40(4), 2021 (additive CCD); closest-point routines
  after C. Ericson, *Real-Time Collision Detection*, 2004.

## Licence

MIT, like the rest of the workspace.

## Comparison with the authors' program

The authors' program (`repulsive-curves`, MIT) is built outside this repository with a headless
batch mode and a fixture driver (`RepulsiveCurves-build/batch-mode.patch`,
`dump-fixtures.patch`; `SPEC.md` §17c). Our energy is the paper's edge-pair discretization and the
reference's is a vertex-pair one (D1), so **energies are never matched**. What is comparable is
compared.

**Matched exactly** (`tests/reference_fixtures.rs`, fixtures in `tests/fixtures/reference/`; trefoil
with 40 vertices, a 99-vertex subsample of the Freedman curve, the full 592-vertex Freedman 600
(digest), an open curve with pins, a junction network, small loops; exponents (3, 6) and (2, 4.5)):

| quantity | agreement |
|---|---|
| dense Gram matrix `A = B + B⁰`, and `B`, `B⁰` separately (midpoint rule, edge tangents) | relative to the largest entry 4e-16 to 5e-15 (tolerance 1e-12); the trapezoid rule differs by 20 % to 200 % |
| edge-length, total-length and point constraint Jacobians | 5e-16 (tolerance 1e-13); our rows are the negatives of the reference's (we differentiate `ℓ⁰ − ℓ`) |
| one coarsening level: surviving vertices, neighbours of the removed ones, prolongation weights, coarse edges (as sets of fine indices; pins as `Point` and `Tangent(Vertex)`) | identical |

**Not comparable**, with reasons in the test file's header: energies and differentials (D1);
barycenter Jacobian (D6; the reference uses constant dual-length weights, checked against its
formula); tangent and surface constraints and penalties (D7, D8, D19, D20); coarsening of closed
loops with fewer than 6 vertices (D17; the test asserts the documented difference), coarsening beyond
level 1 and coarse positions (D17); block-cluster-tree and BVH statistics (D15, D16).

**Compared by final shape** (`examples/compare_reference.rs`, `runs_tmp/rc-wp15/`): the reference's
batch mode (`--no-bh`: exact energy, dense solver; Armijo with the Sobolev–L² start; stop when the
Sobolev dot is below 1e-4) against `flow` with the exact energy and the dense solver, constraints
edge lengths + barycenter, absolute projection tolerance 1e-4, `threads = 5`
(the reference with `OMP_NUM_THREADS = 5`), the guarded line search of the time (commit
`6af096e` plus working tree; the later multigrid, Barnes–Hut adjoint and line-search-default
changes are not in these runs), stop `‖g̃‖₂ < 1e-4`, 5 h wall-time cap. All 12 jobs ran at the
same time on a machine that other jobs were also using (load about 100 to 200), so **times are
approximate**. Curves: Freedman 600 with (2, 4.5) and the first five `perturbed` Knot128 knots
with (3, 6) (about 2000 vertices; dense LU of 6000 unknowns).
Energy: our Eq. 17 evaluated on the output curve. The reference's output is the last OBJ it wrote
(every 10th step).

| curve | vertices | input energy | reference: steps, time, our energy | ours: steps, time, our energy | RMS distance of the curves after arc-length resampling and Procrustes |
|---|---|---|---|---|---|
| Freedman 600 | 592 | 2281 | 504 (converged), 2493 s, 6.0821 | 647 (converged), 2159 s, 6.0815 | 0.0045 (curve length 21.5, bounding box 2.7) |
| Knot128 perturbed 0001 | 1998 | 1652 | 446, cap, 34.54 | 785, cap, 3.258 | not converged |
| 0002 | 1611 | 1638 | 721, cap, 21.86 | 1485, cap (plateau, 1.1935), 1.1935 | not converged |
| 0003 | 1653 | 1569 | 688, cap, 4.324 | 1439, cap (plateau), 1.9958 | not converged |
| 0004 | 2175 | 2315 | 359, cap, 52.54 | 805, cap, 9.906 | not converged |
| 0005 | 2160 | 1982 | 358, cap, 48.96 | 760, cap (plateau), 2.2565 | not converged |

On Freedman 600 both programs reach the same shape and the same energy of our discretization to
1e-4 relative; the remaining distance is a small fraction of the curve size. The five knots did
not converge within 5 hours in the reference, which makes about 1 step/min at this size under
load; ours takes about 2.5 steps/min and larger steps (on 0002, 0003 and 0005 it reaches an energy
plateau with gradient norm below 1e-2 and no decrease in the 8th digit, but the stop threshold
`1e-4` is not met), so the final curves are not the same shape (the knots expand with the fixed edge
lengths, bounding box 130 for the reference and 350 for ours on 0001) and a Procrustes distance is
not meaningful; only energies, steps and times are reported. The loosest honest statement is that
the reference needs more than 2 to 5 times as many hours as the cap to finish these knots.
