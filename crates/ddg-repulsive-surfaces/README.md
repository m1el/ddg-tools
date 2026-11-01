# ddg-repulsive-surfaces

A Rust implementation of

> Chris Yu, Caleb Brakensiek, Henrik Schumacher, Keenan Crane. *Repulsive Surfaces.* ACM
> Transactions on Graphics 40(6), Article 268, 2021. <https://doi.org/10.1145/3478513.3480521>

the tangent-point energy of triangle meshes in R³ (orientation not required: Klein bottles and
Möbius bands are fine), its descent in the fractional Sobolev metric of the paper, constraints
and penalties, Barnes–Hut and hierarchical-matrix acceleration with GMRES, and dynamic remeshing.

**Status: in progress.** The contract is [`SPEC.md`](SPEC.md): every discrete formula with its
`rs-` formula id, where the paper or the reference code is wrong or silent and what we do
(decisions D1–D31), the API, the verification and measurement plans and the work packages
(§15).

Done so far (WP1): the validated mesh type `Surface` (edges, fans, components, boundary loops,
orientation classes), per-face geometry and its derivatives, the cotan stiffness and lumped mass
matrices, the exponents, the tangent-point kernel, the exact energy and its exact differential
(parallel, bitwise independent of the thread count), OBJ input and output.

## Credits and licence

The authors' reference implementation (`repulsive-surfaces`, MIT, "Copyright (c) 2019 Nick
Sharp") may be read and ported; its licence notice is kept in
[`LICENSES/repulsive-surfaces-MIT.txt`](LICENSES/repulsive-surfaces-MIT.txt). The generic
machinery (bounding volume hierarchy, block cluster tree, L-BFGS, implicit surfaces, CCD
primitives) comes from the sibling crate `ddg-repulsive-curves` (Yu, Schumacher, Crane,
*Repulsive Curves*, ACM TOG 40(2), 2021). This crate is MIT licensed, as the workspace.
