# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "safetensors>=0.4"]
# ///
"""Convert the Points as Tori pretrained weights (a pickled Flax NNX State) to safetensors.

The source pickle is the authors' `FundamentalFormPredictor.pkl` (points-as-tori, MIT,
"Copyright (c) 2026 Nicole Feng"); the converted file keeps that notice (see
`crates/ddg-points-as-tori/LICENSES/points-as-tori-MIT.txt`).

Safety: the pickle is read with a restricted `Unpickler` whose `find_class` admits exactly the
seven globals the file uses (SPEC.md §6.2, NOTES-REFERENCE.md §2.4); every other global is an
error. Plain `pickle.load` is never used. Run in isolated mode:

    uv run --no-project --with 'numpy>=2' --with 'safetensors>=0.4' -- \
        python3 -I tools/convert_weights.py <in.pkl> <out.safetensors> [--reference-commit HASH]

Writes `<out.safetensors>` (f32 tensors under their Flax paths, Flax shapes, metadata with the
config, the source sha256, the reference commit and the licence notice) and
`<out.safetensors>.sha256` (sha256 of the output file and of every tensor's bytes).
"""
import argparse
import hashlib
import json
import pickle
import sys

import numpy as np
import numpy._core.multiarray as _ma
from safetensors.numpy import save_file


class _Stub:
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.state = None

    def __setstate__(self, state):
        self.state = state


class State(_Stub):
    """Stand-in for flax.nnx.statelib.State."""


class Param(_Stub):
    """Stand-in for flax.nnx.variablelib.Param."""


class TraceState(_Stub):
    """Stand-in for flax.nnx.tracers.TraceState (discarded)."""


def _reconstruct_array(fun, args, arr_state, aval_state):
    """Stand-in for jax._src.array._reconstruct_array: rebuild the numpy array only."""
    if fun is not _ma._reconstruct:
        raise pickle.UnpicklingError("_reconstruct_array called with an unexpected function")
    arr = fun(*args)
    arr.__setstate__(arr_state)
    return arr


ALLOWED = {
    ("flax.nnx.statelib", "State"): State,
    ("flax.nnx.variablelib", "Param"): Param,
    ("flax.nnx.tracers", "TraceState"): TraceState,
    ("jax._src.array", "_reconstruct_array"): _reconstruct_array,
    ("numpy._core.multiarray", "_reconstruct"): _ma._reconstruct,
    ("numpy", "ndarray"): np.ndarray,
    ("numpy", "dtype"): np.dtype,
}


class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            return ALLOWED[(module, name)]
        except KeyError:
            raise pickle.UnpicklingError(f"forbidden global {module}.{name}") from None


EXPECTED_CONFIG = dict(embed_dim=128, n_freq=0, n_heads=8, n_layers=8, mlp_dim=512, k_neighbors=64, dropout_rate=0.0)


def expected_shapes(cfg):
    e, h, m, L = cfg["embed_dim"], cfg["n_heads"], cfg["mlp_dim"], cfg["n_layers"]
    d = e // h
    s = {
        "input_encoder.mlp.layers.0.kernel": (6, e),
        "input_encoder.mlp.layers.0.bias": (e,),
        "input_encoder.mlp.layers.1.scale": (e,),
        "input_encoder.mlp.layers.1.bias": (e,),
        "final_norm.scale": (e,),
        "final_norm.bias": (e,),
        "mlp.layers.0.kernel": (e, 6),
        "mlp.layers.0.bias": (6,),
    }
    for i in range(L):
        p = f"blocks.{i}"
        for nrm in ("norm1", "norm2"):
            s[f"{p}.{nrm}.scale"] = (e,)
            s[f"{p}.{nrm}.bias"] = (e,)
        for q in ("query", "key", "value"):
            s[f"{p}.attention.{q}.kernel"] = (e, h, d)
            s[f"{p}.attention.{q}.bias"] = (h, d)
        s[f"{p}.attention.out.kernel"] = (h, d, e)
        s[f"{p}.attention.out.bias"] = (e,)
        s[f"{p}.ffn.layers.0.kernel"] = (e, m)
        s[f"{p}.ffn.layers.0.bias"] = (m,)
        s[f"{p}.ffn.layers.3.kernel"] = (m, e)
        s[f"{p}.ffn.layers.3.bias"] = (e,)
    return s


def walk(obj, path, out):
    if isinstance(obj, State):
        inner = obj.state if obj.state is not None else (obj.args[0] if obj.args else {})
        walk(inner, path, out)
        return
    if isinstance(obj, Param):
        st = obj.state
        raw = st.get("_raw_value") if isinstance(st, dict) else st
        if not isinstance(raw, np.ndarray):
            raise ValueError(f"parameter {'.'.join(path)} is not an array ({type(raw).__name__})")
        out[".".join(path)] = raw
        return
    if isinstance(obj, dict):
        for key, val in obj.items():
            if key == "_mapping":
                walk(val, path, out)
            else:
                walk(val, path + [str(key)], out)
        return
    raise ValueError(f"unexpected object at {'.'.join(path)}: {type(obj).__name__}")


NOTICE = (
    "Converted from the pretrained weights of points-as-tori (MIT License, Copyright (c) 2026 Nicole Feng). "
    "The MIT notice in crates/ddg-points-as-tori/LICENSES/points-as-tori-MIT.txt covers this file."
)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pkl")
    ap.add_argument("out")
    ap.add_argument("--reference-commit", default="253bf9fc5e79831a03f0486e0c42feeb67c1c546")
    ap.add_argument("--source", default="points-as-tori src/pointsastori/models/FundamentalFormPredictor.pkl")
    args = ap.parse_args()

    raw_bytes = open(args.pkl, "rb").read()
    src_sha = hashlib.sha256(raw_bytes).hexdigest()
    import io

    data = RestrictedUnpickler(io.BytesIO(raw_bytes)).load()
    if not isinstance(data, dict) or "params" not in data:
        sys.exit("error: the pickle is not a dict with 'params'")
    cfg = {k: data.get(k) for k in EXPECTED_CONFIG}
    if cfg["n_freq"] is None:
        cfg["n_freq"] = 0
    for k, v in EXPECTED_CONFIG.items():
        if cfg[k] != v:
            sys.exit(f"error: config {k} = {cfg[k]!r}, expected {v!r}")

    tensors = {}
    walk(data["params"], [], tensors)
    want = expected_shapes(cfg)
    missing = sorted(set(want) - set(tensors))
    extra = sorted(set(tensors) - set(want))
    if missing or extra:
        sys.exit(f"error: missing tensors {missing}, unexpected tensors {extra}")
    for name, shape in want.items():
        a = tensors[name]
        if a.shape != shape:
            sys.exit(f"error: {name} has shape {a.shape}, expected {shape}")
        if a.dtype != np.float32:
            sys.exit(f"error: {name} has dtype {a.dtype}, expected float32")
        if not np.isfinite(a).all():
            sys.exit(f"error: {name} has non-finite values")
    if len(tensors) != 136:
        sys.exit(f"error: {len(tensors)} tensors, expected 136")

    out = {k: np.ascontiguousarray(tensors[k]) for k in sorted(tensors)}
    n_params = int(sum(a.size for a in out.values()))
    meta = {
        "config": json.dumps({k: cfg[k] for k in ("embed_dim", "n_layers", "n_heads", "mlp_dim", "k_neighbors", "n_freq")}),
        "source": args.source,
        "source_sha256": src_sha,
        "reference_commit": args.reference_commit,
        "license": NOTICE,
        "format": "flax-nnx paths, flax shapes, float32 little-endian",
    }
    save_file(out, args.out, metadata=meta)
    file_sha = hashlib.sha256(open(args.out, "rb").read()).hexdigest()
    with open(args.out + ".sha256", "w") as f:
        f.write(f"{file_sha}  {args.out.rsplit('/', 1)[-1]}\n")
        f.write(f"# source {args.source} sha256 {src_sha}\n")
        for k, a in out.items():
            f.write(f"# tensor {k} {list(a.shape)} {hashlib.sha256(a.tobytes()).hexdigest()}\n")
    print(f"wrote {args.out}: {len(out)} tensors, {n_params} parameters, sha256 {file_sha}")


if __name__ == "__main__":
    main()
