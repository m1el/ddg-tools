"""WP10: criterion output (the two bench logs) -> markdown table (median time per benchmark)."""
import re, sys, collections
def parse(path):
    out = collections.OrderedDict(); name = None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"^(representation|retriangulation)/(.*)$", line)
        if m and "Benchmarking" not in line: name = line.strip(); continue
        m3 = re.search(r"time:\s+\[([\d.]+) (\S+) ([\d.]+) (\S+) ([\d.]+) (\S+)\]", line)
        if m3:
            key = line[:m3.start()].strip() or name
            if key: out[key] = (float(m3.group(3)), m3.group(4))
    return out
def ms(v): return v[0] * {"ns": 1e-6, "µs": 1e-3, "ms": 1.0, "s": 1e3}[v[1]]
r = parse(sys.argv[1]); s = parse(sys.argv[2]) if len(sys.argv) > 2 else {}
print("| benchmark | median |"); print("|---|---|")
for k, v in r.items(): print(f"| {k} | {ms(v):.3f} ms |")
print("\nretriangulation (ddg-intrinsic):\n"); print("| benchmark | median |"); print("|---|---|")
for k, v in s.items(): print(f"| {k} | {ms(v):.3f} ms |")
