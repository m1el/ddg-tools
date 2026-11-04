# ddg-points-as-tori

Pointwise signed distance for oriented point clouds, after

> Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed Distance
> for Point Clouds.* ACM Transactions on Graphics 45(4), Article 53, 2026.
> <https://doi.org/10.1145/3811385>

A small transformer predicts, for each point, the coefficients of a local quadratic surface; each
point gets a torus with the same shift, principal curvatures and directions, and the signed
distance at any query is a distance-weighted blend of the tori's closed-form SDFs (Eq. 23–26).

Status: contract only (`SPEC.md`); the code is written by the work packages listed there.
`NOTES-REFERENCE.md` is the survey of the authors' reference implementation.

## Credits and licences

- The method and the pretrained network are by Nicole Feng, Ioannis Gkioulekas and Keenan Crane.
- The reference implementation `points-as-tori` (MIT, "Copyright (c) 2026 Nicole Feng") is read
  and ported where noted; its notice is in `LICENSES/points-as-tori-MIT.txt`. The converted
  pretrained weights shipped with this crate come from that repository and are under the same
  notice.
- This crate: MIT (workspace licence).
