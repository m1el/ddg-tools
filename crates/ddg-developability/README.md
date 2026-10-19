# ddg-developability

An implementation, from the text of the paper (clean room), of

> Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.*
> ACM Transactions on Graphics 37(4), Article 77, 2018. <https://doi.org/10.1145/3197517.3201303>

A vertex of a triangle mesh is *discrete developable* (a hinge) when its triangles split into two
flat, edge-connected pieces — equivalently, when its triangle normals lie in one plane. The crate
evaluates the paper's six per-vertex measures of how far a vertex is from a hinge (combinatorial
width and its max variant, normal covariance and its max variant, the intrinsic width in the log
map of the Gauss sphere and its branching-free variant), their gradients, and drives a mesh
toward a piecewise developable surface with L-BFGS and the Lewis–Overton line search, 4-1
subdivision between rounds and small-angle remeshing. A cut through the remaining seam vertices
and a Boundary First Flattening (`ddg-bff`) turn the result into near-isometric UV charts.

The contract — every formula with its derivation, the differences from the paper as printed
(`D1`–`D29`), the API, the formula ids and the verification plan — is [`SPEC.md`](SPEC.md).

The authors' implementation is GPL and is not used, read or paraphrased here; the 3x3
eigensolver is a textbook cyclic Jacobi method, not Kopp's LGPL code. MIT licensed, like the
rest of the workspace.
