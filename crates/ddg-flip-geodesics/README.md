# ddg-flip-geodesics

Geodesic paths, loops and networks on triangle meshes by intrinsic edge flips, after

> Nicholas Sharp, Keenan Crane. *You Can Find Geodesic Paths in Triangle Meshes by Just
> Flipping Edges.* ACM Transactions on Graphics 39(6), Article 249, 2020.
> <https://doi.org/10.1145/3414685.3417839>

run on either intrinsic triangulation of this workspace: the signposts of `ddg-intrinsic`
(Sharp, Soliman, Crane 2019) and the integer coordinates of `ddg-integer-coords` (Gillespie,
Sharp, Crane 2021). Everything is generic over `ddg_intrinsic::Triangulation`.

The contract of the crate is [`SPEC.md`](SPEC.md): scope, data, the geometry of a joint, the
algorithm, verification, and every decision the paper leaves open (`F1`–`F25`).

## What is implemented (work package F1)

- Edge networks: paths and loops as sequences of directed edges (`Segment`, `EdgePath`),
  `closed` explicit, pinned vertices (path ends pinned implicitly), multiply covered edges with
  the layering of the segments that share an edge (paper Fig. 11; default as the reference,
  or given explicitly), stable segment ids.
- Validity (`Network::check`): connected paths, closed loops, pinned ends, a complete layering
  and non-crossing at every vertex (the chords of the joints, in the cyclic order of the
  segment ends, do not interleave). Crossing inputs are refused, never repaired.
- Initial paths: through a vertex sequence, and Dijkstra on the edge graph of the current
  triangulation with ties by `(distance, halfedge id)`.
- The geometry of a joint: both wedge angles summed corner by corner (Eq. 1), boundary wedges
  (`∞`), the Δ-complex cases (a self-edge incident on a wedge twice, the out-and-back joint at a
  vertex of degree one with a wedge of degree zero), the classification straight / left /
  right with the reference's tolerance, flexible wedges and wedges blocked by a path end.
- The geodesic certificate on the triangulation (SPEC §7 layer 1a).
- The options and report types of the shortening.

Next: FlipOut (Alg. 1, 2, loops, networks), extraction on the input mesh for both
representations, geodesic Bézier curves, the intrinsic constrained Delaunay triangulation
(SPEC §12).

## Credits

The algorithm is the paper's. The authors' implementation in geometry-central (MIT,
© 2017–2019 Nicholas Sharp and the geometry-central contributors; notice in
[`LICENSES/geometry-central-MIT.txt`](LICENSES/geometry-central-MIT.txt)) was read for the
details the paper leaves out; the places that follow it are cited `gc: <file>:<line>` in the
docs and in `SPEC.md`.
