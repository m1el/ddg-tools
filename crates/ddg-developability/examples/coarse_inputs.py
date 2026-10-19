# /// script
# dependencies = ["pymeshlab", "numpy"]
# ///
"""Coarse, uniform inputs for the paper's workflow (QUALITY.md; paper Sec. 5.2: the inputs of Figs. 1,
22-25 were remeshed with Instant Meshes to uniform triangles, about 1k for Fig. 1).

We have no Instant Meshes here; MeshLab's isotropic explicit remeshing (pymeshlab, a tool used only to
make test inputs, not part of the crate) gives uniform, nearly equilateral triangles of a chosen edge
length on the dense model.

    uv run crates/ddg-developability/examples/coarse_inputs.py OUTDIR [TARGET_FACES ...]

Writes OUTDIR/<model>_iso<N>.obj for spot, bob, blub (dense triangulated models of data/models).
"""

import math
import os
import sys

import numpy as np
import pymeshlab as ml

MODELS = {
    "spot": "data/models/spot/spot/spot_triangulated.obj",
    "bob": "data/models/bob/bob/bob_tri.obj",
    "blub": "data/models/blub/blub/blub_triangulated.obj",
}


def area(m):
    V, F = m.vertex_matrix(), m.face_matrix()
    return 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum()


def main():
    out = sys.argv[1]
    targets = [int(x) for x in sys.argv[2:]] or [1000]
    os.makedirs(out, exist_ok=True)
    for name, path in MODELS.items():
        for n in targets:
            ms = ml.MeshSet()
            ms.load_new_mesh(path)
            A = area(ms.current_mesh())
            L = math.sqrt(4 * A / (math.sqrt(3) * n))
            ms.meshing_isotropic_explicit_remeshing(targetlen=ml.PureValue(L), iterations=15)
            m = ms.current_mesh()
            dst = os.path.join(out, f"{name}_iso{n}.obj")
            ms.save_current_mesh(dst, save_vertex_normal=False, save_textures=False, save_wedge_texcoord=False,
                                 save_vertex_color=False)
            print(f"{dst}: {m.vertex_number()} vertices, {m.face_number()} faces, target edge {L:.4f}")


if __name__ == "__main__":
    main()
