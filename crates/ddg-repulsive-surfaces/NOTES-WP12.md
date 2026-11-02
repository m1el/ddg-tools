# WP12 notes (generators, scene reader, fetch script)

Build/test: `CARGO_TARGET_DIR=$PWD/target-wp12 cargo test --release -p ddg-repulsive-surfaces --test generate --test scene`
(about 6 s). Data: `scripts/rs-fetch-scenes.sh` copies the reference scenes to `data/repulsive-surfaces/scenes/`.

## Scene reader (`src/scene.rs`)
- `SceneFile::parse/read` (syntax, no I/O, `to_text()` round trip) and `SceneSpec::resolve/load` (mesh, obstacles read;
  penalties built). Errors are `RsError::Parse { path, line, message }` (line 0: whole file); a missing file is `RsError::Io`
  naming it.
- Stricter than the reference: unknown keywords (incl. the `minimzie` typo), surplus tokens, missing `α β` (the reference leaves
  them uninitialized: `Boundary/scene-tall.txt` is the one reference scene that does not parse), `constrain area 3` without
  iterations, duplicate `repel_mesh/method/log/iteration_limit/time_limit/autotarget_volume`. Whitespace-split (tabs ok).
  `autotarget_volume` without a number means 1. `iteration_limit 0` / `time_limit 0` = no limit (`None`).
- Mapping: `minimize area|volume w` = `TotalArea|TotalVolume { power: 2 }`; `area_deviation|volume_deviation` targets = the
  initial area / signed volume; `squared_error` targets = start positions; `boundary_length w L` = `BoundaryLength{target: L}`
  (paper's squared form, D24); `boundary_curvature` = `dual_length: false`; `willmore`; obstacles = `MeshObstacle`
  (`ObstacleKind::TangentPoint`, leaf 8) / `PointCloudObstacle`, `centered` subtracts the area-weighted barycenter; implicit
  = `ImplicitObstacle|ImplicitAttractor` with `Quadrature::PerVertex`; torus = `FramedTorus` with axis `y`; attractor set =
  vertices with a corner UV `u > 1e-10` (None = all when the OBJ has no such corner).
- `SceneSpec` fields for WP5: `constraints: Vec<SceneConstraint { constraint: Constraint, schedule: ScheduleSpec }>`
  (`Growth {multiplier, iterations, addition}` -> `Schedule`; `Move {displacement, iterations}` -> `PinSchedule`, every pin of
  the statement gets the displacement); `normal_pins` (`constrain boundary_normals|normals`: the library has no normal
  constraint, WP5 reports `Unsupported` or implements it); `penalties` in the reference's order (potentials, obstacles,
  implicits) with `penalty_lines`; `obstacle_volume` + `autotarget_volume` (final volume target = ratio * obstacle_volume);
  `method: SceneMethod` (`willmore` = replace the energy, H2 flow; `l2`, `bqn` = BCQN ...), `allow_barycenter_shift`,
  `disable_near_field` (= `low_order: false`), limits, `log`. `exponents`, `surface`, `positions`.
  Each `constrain vertices` line is its own `Pins` constraint (the reference merges them).
- Reference scenes: 36 scene files; 25 resolve fully; 10 refer to missing meshes (`figure-eight-smoothed.obj` x2, `hshape.obj`
  x3, `colliding.obj`, `remeshed.obj`, `doubletorus4.6k.obj`, `SlideOffBar5.obj`, `end.obj`, `boxscene.txt`: `box.obj` and
  `doubletorus1.6k.obj`); they parse; 1 (`scene-tall.txt`) is an error as above.

## Generators (`src/generate/`)
- `tetra`: marching tetrahedra (Kuhn 6-tet cubes, vertices on grid edges, outward orientation from the field sign), block
  skipping with the Lipschitz bound, vertices clamped >= 0.15 cell from a grid vertex, 6 rounds of Laplacian smoothing each
  followed by two Newton steps onto the zero set. `mesh_with_faces(field, box, target)` adjusts the cell width (within 1-2 %
  of the target in practice). Not the WP4 remesher (not merged): callers remesh if they need.
- `shapes`: `genus_surface(g, faces)` (rounded slab with g holes in a row, g = 0..11 tested); `handcuffs(k, faces)`: tube
  (genus 2) around the cores of the torus link T(2,2k) (linking number k, checked by the Gauss integral) plus a bridge;
  `knotted_tunnel(Knot, faces)`: a flat-topped ball minus a tube along a trefoil / figure-eight knot, the lobe sticking out of
  the top chosen by searching 400 directions (genus 1, one component, a knotted torus); `fig12(i)`: the five D31 analogues:
  25 076 / 26 348 / 21 640 / 44 848 / 23 344 faces (paper 25 098 / 26 306 / 21 712 / 45 198 / 23 206). The tube is
  resolved only above about 20 000 faces for k = 6 and 30 000 for the figure-8 tunnel; coarser requests return an error
  (topology check), not a wrong surface. `graph_tube(curves, radius, faces)` meshes any tube network.
  Interpretation of "tunnel": the paper's text was not consulted for the exact figure; this is our reading, an analogue.
- `primitives`: sphere, torus, box, capsule, punctured torus (genus 1, one boundary loop, chi = -1), `cobblestones`
  (n seeded spheres inside a slab obstacle), `soup_with_holes(mesh, fraction, seed)`, `h_shape_coarse()` (32 vertices, 60
  faces; vertices 0..23 are the prism corners = the 24 pins of the `CoarseInterpolation` scenes, 24..31 free).
- `perturbed::PerturbedTorus`: `f0 + Phi(f0)`, `Phi` a seeded trigonometric polynomial R^3 -> R^3 (our choice; the paper does not
  give it), `sample(nu, nv, squeeze_u, squeeze_v)` with a sine warp of the parameter grid.
- `check`: brute-force sweep self-intersection test, Gauss linking number, clearance.

## Not done
- Ishii et al. 2012 Table 1 (Fig. 16): the paper is not available offline and its diagrams cannot be transcribed reliably
  from memory; nothing was invented. 4_1 and 5_x are missing. `graph_tube` is the entry point once diagrams are transcribed
  (polylines of the trivalent graph).
- Hopf-linked handcuffs (k = 1) are a *trivial* handlebody-knot (its complement is a genus-2 handlebody); they are not 4_1.

## SPEC proposals
- Section 14.4: say that the Fig. 16 table is blocked on access to Ishii et al.; Fig. 22's coarse H has 24 pinned + 8 free
  vertices in our analogue; D22 note that the loader computes the attractor set (`vertices: Some(..)`).
- Section 9.1: `Parse { path: Option<String>, line, message }` is the implemented form.
