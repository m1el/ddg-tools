# Repulsive Surfaces: applications of Sec. 9

Rust example `apps` (`crates/ddg-repulsive-surfaces/examples/apps/`) runs every application; the
smoke versions are `crates/ddg-repulsive-surfaces/tests/applications.rs`.

```
scripts/rs-fetch-scenes.sh                       # reference scenes -> data/repulsive-surfaces/scenes/
showcase/repulsive-surfaces/run_all.sh           # all applications, 16 jobs x 4 threads, OBJ frames
showcase/repulsive-surfaces/run_all.sh --smoke   # reduced sizes, no output files
target/release/examples/apps list                # names
target/release/examples/apps genus05 --out DIR --frame-every 10 --threads 8
```

Frames go to `showcase/media/repulsive-surfaces/<name>/frame_NNNNN.obj` (untracked), the per-step
CSV to `<name>/steps.csv`, logs to `runs_tmp/wp15/<mode>/<name>.log` (`RESULT` lines: steps,
time, energy, area/volume, constraint violation, self-intersection by `ccd::intersects`).
The eversion writes `eversion/evert_NNNNN.obj`: the forward flow, then the reflected and
rotated frames in reverse order. Results and open problems: `crates/ddg-repulsive-surfaces/NOTES-WP15.md`.

| name | paper | what |
|---|---|---|
| `genusNN` | Fig. 15 | genus 0..11 from a linear arrangement, fixed area and barycenter; moment-tensor symmetry indicators |
| `genusNN_plane5`, `_plane100` | Fig. 17 | with an attracting plane (weight 5; 100 as the strong stand-in for a plane constraint) |
| `handcuffs` | Fig. 1, 14 | linked and unlinked handcuffs, plane attractor, same area; final energies compared |
| `pole`, `slideoff` | Fig. 18, 19 | handle off a pole / bar: implicit cylinder / torus, plane attractor, free barycenter; stage 2 with the repulsion weakened |
| `eversion` | Fig. 20 | punctured torus, zero-volume penalty, fixed area, boundary length and curvature penalties |
| `h_*` | Fig. 21, 22 | H with sparse/dense pins, hard pins or squared error, volume/area growth, with and without a box |
| `shrinkwrap_*` | Fig. 23 | sphere shrunk onto a point cloud (bunny, spot) and onto a polygon soup with 35 % of the faces deleted |
| `envelopes` | Fig. 24 | nested envelopes, volume x1.15 per level, previous level as obstacle, coarser each level |
| `walnut`, `cobblestones` | Fig. 25 | area growth inside a shell; 14 spheres growing inside a slab |
| `subcritical`, `klein` | Fig. 26 | exponents below the critical ones with `low_order: false` |
