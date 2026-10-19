# /// script
# dependencies = []
# ///
"""Markdown table of a quality study directory (examples/quality_study.sh): one row per job log."""
import glob
import os
import sys

COLS = ["faces_out", "e_final", "status", "gnorm_final", "round_iterations", "seconds", "bad_faces",
        "e30", "e30_frac", "folds", "seam_comp", "seam_short", "seam_mean_len", "seam_long_frac",
        "d_mean", "d_haus", "vol_ratio", "min_angle", "faces_lt1"]


def parse(path):
    d = {}
    for line in open(path, errors="replace"):
        if line.startswith("RESULT ") or line.startswith("QUALITY "):
            for kv in line.split()[1:]:
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = v
    return d


def main():
    out = sys.argv[1]
    cols = sys.argv[2].split(",") if len(sys.argv) > 2 else COLS
    rows = []
    for p in sorted(glob.glob(os.path.join(out, "*.log"))):
        d = parse(p)
        rows.append([os.path.basename(p)[:-4]] + [d.get(c, d.get("error", "")[:30] if c == "status" else "") for c in cols])
    print("| tag | " + " | ".join(cols) + " |")
    print("|---" * (len(cols) + 1) + "|")
    for r in rows:
        print("| " + " | ".join(r) + " |")


if __name__ == "__main__":
    main()
