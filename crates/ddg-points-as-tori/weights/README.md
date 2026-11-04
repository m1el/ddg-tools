# Pretrained weights

`pat-reference.safetensors` is the pretrained network of

> Nicole Feng, Ioannis Gkioulekas, Keenan Crane. *Points as Tori: Fast Pointwise Signed Distance
> for Point Clouds.* ACM Transactions on Graphics 45(4), Article 53, 2026.

converted from the authors' reference implementation, `points-as-tori` at commit
`253bf9fc5e79831a03f0486e0c42feeb67c1c546`, file `src/pointsastori/models/FundamentalFormPredictor.pkl`
(sha256 `aaa2f3882c36ac32d57361377a73ac830c47d03d6167d51249127a6c40ef7ff9`).

**Licence.** The weights are part of `points-as-tori`, MIT License, Copyright (c) 2026 Nicole Feng.
The notice is in `../LICENSES/points-as-tori-MIT.txt` and covers this file; the same notice is in
the file's safetensors metadata (`license`).

**Format.** safetensors: 136 f32 tensors, 1,588,358 parameters, named by their Flax NNX paths
(`blocks.3.attention.query.kernel`, …) with Flax shapes (`[128, 8, 16]` attention kernels, `[in, out]`
linear kernels). The metadata holds `config` (embed 128, 8 layers, 8 heads, MLP 512, k 64, no
frequencies), `source`, `source_sha256` and `reference_commit`. `pat-reference.safetensors.sha256`
lists the file hash and a hash per tensor.

**How it was made.** `tools/convert_weights.py`, which reads the pickle with a restricted unpickler
that admits only the seven globals the file uses (no plain `pickle.load`), checks the configuration,
the 136 names, shapes, dtypes and finiteness, and writes the file:

```
uv run --no-project --with 'numpy>=2' --with 'safetensors>=0.4' -- \
    python3 -I tools/convert_weights.py \
    /var/data/discrete-geometry/PointsAsTori/points-as-tori/src/pointsastori/models/FundamentalFormPredictor.pkl \
    weights/pat-reference.safetensors
```

The converted tensors were compared with the parameters the authors' own loader produces
(`runs_tmp/pat-ref/activations_tooth/weights.npz`, written by the reference build): all 136 equal,
max |Δ| = 0.

The crate embeds the file with `include_bytes!` behind the default feature `bundled-weights`
(`Model::reference()`).

The second pickle of the training data (`FundamentalFormPredictor_pretrained.pkl`, a different set of
weights with the same configuration, NOTES-REFERENCE.md §3) is converted with the same tool into
`data/points-as-tori/` and is not committed.
