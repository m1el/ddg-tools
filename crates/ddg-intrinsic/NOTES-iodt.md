# Working notes: vertex repositioning and iODT (checkpoint; fold into SPEC at the end)

Branch: worktree of `intrinsic-triangulations` (fast-forwarded to a0773fd). Build dir `target-iodt`.
`data/` holds symlinks to the main checkout's data.

## Decisions (to become SPEC N19..)

* `move_vertex(v, [dx, dy])`: displacement in the vertex's own frame (signpost frame). Inserted
  vertices only. Validity: exact kernel test (every new triangle positively oriented with the
  shape limit of N5), via `2A' = 2A - cross(v, j_{n+1} - j_n)` in the star layout; the paper's
  height test is sufficient, not necessary. Outside: error with the fraction at which a new
  triangle degenerates. Never shortened by `move_vertex`.
* Lengths: `l'^2 = (l - r cos a)^2 + (r sin a)^2` = law of cosines (Alg. 10).
* Place: trace over the input mesh from the vertex's old place along v (Alg. 10 traces from a
  neighbour j0; same point in exact arithmetic, shorter path here). Snapping as insertion.
* Signposts: outgoing anchored on the trace's back direction + layout angle, then cumulative
  corner angles (vertex update); incoming by signpost update from the link halfedge.
* `edge_on_input` of the incident edges cleared on an interior move.
* Boundary inserted vertex: moves along the boundary only (frame x-axis; v[1] must be 0),
  place by interpolation on its input edge.
* Self-edge at the vertex: refused.
* iODT large step: Chen & Holst line search (alpha = step, halved until accepted, at most
  `max_halvings`), counted. Paper's general route (remove/retriangulate/reinsert) not done.
* iODT centres: circumcentres (paper) or Chen & Holst's barycentres for triangles with a
  boundary vertex (option).

## Cited
* Chen & Holst 2011 preprint downloaded: `../NavigatingIntrinsicTriangulations/cited/chen-holst-2011/ChHo10a.pdf`
  (Lemma 4.2 Eq. 4.5, Eq. 4.13 with line search alpha = 1, 1/2, ...; boundary barycentres Sec. 4).
* Planar identity for the proof (per triangle T = (x, b, c) CCW, rot(u) = (-u_y, u_x)):
  `2x|T| + (|b|^2+|c|^2) grad_x|T| = 2|T|(x - c_T) + F(c) - F(b)`,
  `F(y) = 1/2 (|y|^2 rot(y - x) + |x|^2 rot(y))`; telescopes over a closed star.

## Progress
- [ ] reposition.rs (move_vertex) + tests
- [ ] odt.rs (optimal_delaunay) + tests
- [ ] stats (angle/area summaries)
- [ ] Lean proofs (subagent), oracle fixture (subagent), Rust replay
- [ ] study example, bench
- [ ] docs (SPEC, README, ROADMAP)
- [ ] ddg-uv integration + flow comparison
