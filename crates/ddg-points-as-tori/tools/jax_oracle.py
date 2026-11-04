# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = ["jax==0.8.2", "flax==0.12.2", "numpy>=2", "safetensors>=0.4"]
# ///
"""JAX CPU oracle for the Points as Tori network (SPEC.md §15.1-4).

The architecture of the reference's `FundamentalFormPredictor` (points-as-tori, `network.py`, MIT,
Copyright (c) 2026 Nicole Feng) re-declared with Flax NNX 0.12.2 / JAX 0.8.2 on the CPU, in f32
(`jax_enable_x64 = False`, as the reference). It does not import the reference. The weights come
from our safetensors file (`weights/pat-reference.safetensors`), so the oracle checks the
converted weights too.

Inputs: `points.npy`, `normals.npy` (float64, (n, 3)) and `knn.npy` (int, (n, k), row i starting
with i; e.g. the golden `knn_indices.npy` or the crate's fitting neighborhoods). Output: an
uncompressed `.npz` with `raw` (n, 6) f32 network outputs, `sigma` (n,) f32 and `coeffs` (n, 6) f32
(the reference's rescaling).

    uv run --no-project --python 3.12 --with jax==0.8.2 --with flax==0.12.2 --with 'numpy>=2' \
        --with 'safetensors>=0.4' -- python3 -I tools/jax_oracle.py \
        weights/pat-reference.safetensors points.npy normals.npy knn.npy out.npz [--chunk 4096]
"""
import argparse
import os

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
from flax import nnx  # noqa: E402
from safetensors.numpy import load_file  # noqa: E402

jax.config.update("jax_enable_x64", False)
NAN_EPSILON = 1e-10


def orthonormal_basis(n):
    n_hat = n / jnp.maximum(jnp.linalg.norm(n, axis=-1, keepdims=True), NAN_EPSILON)
    ax = jnp.abs(n_hat[..., 0:1])
    zero, one = jnp.zeros_like(ax), jnp.ones_like(ax)
    u = jnp.where(ax > 0.9, jnp.concatenate([zero, one, zero], -1), jnp.concatenate([one, zero, zero], -1))
    s = jnp.cross(n_hat, u)  # the reference's Rodrigues rotation by 0 is the identity
    s = s / jnp.maximum(jnp.linalg.norm(s, axis=-1, keepdims=True), NAN_EPSILON)
    return s, jnp.cross(n_hat, s)


class Block(nnx.Module):
    def __init__(self, e, h, m, *, rngs):
        self.attention = nnx.MultiHeadAttention(
            num_heads=h, in_features=e, qkv_features=e, out_features=e, decode=False, rngs=rngs
        )
        self.ffn = nnx.Sequential(
            nnx.Linear(e, m, rngs=rngs), nnx.gelu, nnx.Dropout(0.0, rngs=rngs), nnx.Linear(m, e, rngs=rngs),
            nnx.Dropout(0.0, rngs=rngs),
        )
        self.norm1 = nnx.LayerNorm(e, epsilon=1e-6, rngs=rngs)
        self.norm2 = nnx.LayerNorm(e, epsilon=1e-6, rngs=rngs)

    def __call__(self, x):
        x = x + self.attention(self.norm1(x))
        return x + self.ffn(self.norm2(x))


class InputFeatures(nnx.Module):
    def __init__(self, e, *, rngs):
        self.mlp = nnx.Sequential(nnx.Linear(6, e, rngs=rngs), nnx.LayerNorm(e, epsilon=1e-6, rngs=rngs))

    def __call__(self, positions, normals):
        nn_ = normals / jnp.maximum(jnp.linalg.norm(normals, axis=-1, keepdims=True), NAN_EPSILON)
        cn = nn_[:, 0:1, :]
        dv = positions - positions[:, 0:1, :]
        median = jnp.median(jnp.sqrt(jnp.sum(dv * dv, -1)), axis=-1)
        npos = positions / median[:, None, None]
        s_hat, t_hat = orthonormal_basis(cn)
        dv = npos - npos[:, 0:1, :]
        coords = [jnp.sum(dv * s_hat, -1, keepdims=True), jnp.sum(dv * t_hat, -1, keepdims=True), jnp.sum(dv * cn, -1, keepdims=True)]
        local = [jnp.sum(nn_ * s_hat, -1, keepdims=True), jnp.sum(nn_ * t_hat, -1, keepdims=True), jnp.sum(nn_ * cn, -1, keepdims=True)]
        return self.mlp(jnp.concatenate(coords + local, -1)), median


class Net(nnx.Module):
    def __init__(self, e=128, L=8, h=8, m=512, *, rngs):
        self.input_encoder = InputFeatures(e, rngs=rngs)
        self.blocks = nnx.List([Block(e, h, m, rngs=rngs) for _ in range(L)])
        self.final_norm = nnx.LayerNorm(e, epsilon=1e-6, rngs=rngs)
        self.mlp = nnx.Sequential(nnx.Linear(e, 6, rngs=rngs))

    @nnx.jit
    def __call__(self, positions, normals):
        x, sigma = self.input_encoder(positions, normals)
        for b in self.blocks:
            x = b(x)
        raw = self.mlp(self.final_norm(x)[:, 0, :])
        coeffs = jnp.stack([raw[:, 0] * sigma, raw[:, 1], raw[:, 2], raw[:, 3] / sigma, raw[:, 4] / sigma, raw[:, 5] / sigma], -1)
        return raw, sigma, coeffs


def load(model, path):
    w = load_file(path)
    for k, v in w.items():
        obj = model
        parts = k.split(".")
        for p in parts[:-1]:
            obj = obj[int(p)] if p.isdigit() else getattr(obj, p)
        getattr(obj, parts[-1]).value = jnp.asarray(v)
    return len(w)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("weights")
    ap.add_argument("points")
    ap.add_argument("normals")
    ap.add_argument("knn")
    ap.add_argument("out")
    ap.add_argument("--chunk", type=int, default=4096)
    a = ap.parse_args()
    model = Net(rngs=nnx.Rngs(0))
    nt = load(model, a.weights)
    nnx.eval_mode(model)
    P, N, idx = np.load(a.points), np.load(a.normals), np.load(a.knn).astype(np.int64)
    if (idx[:, 0] != np.arange(len(P))).any():
        raise SystemExit("error: knn rows must start with their own index")
    outs = [[], [], []]
    for s in range(0, len(P), a.chunk):
        r = model(jnp.asarray(P[idx[s : s + a.chunk]], jnp.float32), jnp.asarray(N[idx[s : s + a.chunk]], jnp.float32))
        for o, v in zip(outs, r):
            o.append(np.asarray(v))
    raw, sigma, coeffs = (np.concatenate(o) for o in outs)
    np.savez(a.out, raw=raw, sigma=sigma, coeffs=coeffs)
    print(f"{a.out}: {len(P)} neighborhoods of {idx.shape[1]} tokens, {nt} tensors loaded")


if __name__ == "__main__":
    main()
