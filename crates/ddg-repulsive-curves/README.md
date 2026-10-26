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
