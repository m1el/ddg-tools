# Optimizer, subdivision, remeshing: progress notes

Done (SPEC §4, §5.2, §5.3, §10.5, §10.6):

* `src/optimize/`: generic `minimize(objective, x0, Problem{length_scale, group}, opts, callback)`
  (`core.rs`), options/report/statuses (`options.rs`), mesh entry point `optimize` (`mesh.rs`).
  L-BFGS (two-loop, H0 from the newest pair), Lewis-Overton Algorithm 2.6, gradient descent with
  weak Wolfe or backtracking, fixed step, units of the mean edge length, step cap, pinned vertices.
* `src/subdivide.rs`: 4-1 midpoint subdivision, `Subdivided` + `pinned_after`.
* `src/remesh.rs`: flips of caps, collapses of needles, all legality rules, refusal reasons are
  `&'static str` and unit-tested one by one (`cargo test --lib remesh`).
* Tests: `tests/optimize.rs`, `tests/subdivide.rs`, `tests/remesh.rs`.

Decisions that go beyond the SPEC text (all only restrict or clarify):

* Objective trait: `EvalError::Undefined` (+inf in a line search, error at the start) vs `Fatal`.
  The mesh objective maps `UndefinedEnergy`, `DegenerateFace`, `NonFinitePosition` at a trial point
  to undefined.
* Fixed step: an undefined trial point ends the run with `LineSearchFailed`.
* L-BFGS with `Backtracking` is allowed (pair skipped if s^T y <= 0).
* Line-search failure with non-empty memory: reset, retry once from steepest descent; a failure
  from steepest descent ends the run. A non-descent L-BFGS direction resets the memory; a
  non-descent steepest direction ends the run (`NonDescent`) at once.
* Remesh: flips are also refused when a new face's normal turns by >= 90 degrees against an old
  face of the quad (reflex quad), the adjacency rule is tested before the valence rule (the
  valence rules are implied by adjacency / link condition for interior vertices and are kept as
  safety nets), a pinned vertex survives a collapse (SPEC: lower id) and `Remeshed::pinned` gives
  the new pinned ids, `RemeshOptions::validate_each_operation` (tests), `energy_before/after` are
  `Option<f64>` filled in by the flow driver.

Observations:

* BFGS on Nesterov's second Chebyshev-Rosenbrock function (abs inside) stalls at non-minimizing
  Clarke stationary points (known, Gurbuzbalaban-Overton); on the first (squared) variant, n = 3,
  L-BFGS reaches f < 1e-4 in 2000 iterations from 20/20 random starts.

Remains: the flow driver (owner: other agent) must call `remesh` between `optimize` calls
(break from the progress callback, remesh, restart; memory is fresh per call).
