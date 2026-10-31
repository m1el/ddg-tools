"""Cross-tabulates thingi_convert's table with Thingi10k's geometry_data.csv (WP10)."""
import csv, collections, sys
tsv, meta = sys.argv[1], sys.argv[2]
m = {r["file_id"]: r for r in csv.DictReader(open(meta))}
rows = list(csv.DictReader(open(tsv), delimiter="\t"))
c = collections.Counter(); ex = collections.defaultdict(list)
for r in rows:
    g = m.get(r["id"])
    mm = "no metadata" if g is None else ("manifold" if g["vertex_manifold"] == "1" and g["edge_manifold"] == "1" else "not manifold")
    k = (mm, r["status"] if r["status"] == "kept" else "excluded: " + r["reason"])
    c[k] += 1
    if len(ex[k]) < 6: ex[k].append(r["id"])
for k, v in sorted(c.items()): print(f"{v:6d}  metadata {k[0]:13s} ours {k[1]:40s} e.g. {','.join(ex[k])}")
print("metadata manifold:", sum(v for k, v in c.items() if k[0] == "manifold"))
print("kept:", sum(v for k, v in c.items() if k[1] == "kept"))
