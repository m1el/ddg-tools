"""Compose the panels of render/developable_paper_stills.sh (work/dev2/p) into the media figures.
Needs ImageMagick (`convert`). Run with plain python3 from showcase/:
    python3 render/developable_paper_compose.py hero|inputs|energies|crumpled|all
"""
import json
import subprocess
import sys
from pathlib import Path

SHOW = Path(__file__).resolve().parents[1]
W, P, M = SHOW / "work", SHOW / "work/dev2/p", SHOW / "media"
BG, INK, FONT = "#f3f4f6", "#1f2329", "DejaVu-Sans"


def cv(*args):
    subprocess.run(["convert", *map(str, args)], check=True)


def cell(png, out, size, caption=None, ps=None):
    """A panel flattened on the background, resized to `size`, with a caption at the top left."""
    args = [png, "-background", BG, "-alpha", "remove", "-alpha", "off", "-resize", f"{size}x{size}"]
    if caption:
        args += ["-font", FONT, "-pointsize", ps or size // 36, "-fill", INK, "-gravity", "northwest",
                 "-annotate", f"+{size // 40}+{size // 40}", caption]
    cv(*args, out)


def text_block(lines, width, out, ps, bold_first=True):
    args = ["-size", f"{width}x", "-background", BG, "-fill", INK, "-font", FONT, "-pointsize", ps]
    parts = []
    for i, l in enumerate(lines):
        f = out.with_name(out.stem + f"_l{i}.png")
        cv("-background", BG, "-fill", INK, "-font", "DejaVu-Sans-Bold" if (i == 0 and bold_first) else FONT,
           "-pointsize", ps, f"label:{l}", f)
        parts.append(f)
    cv(*parts, "-background", BG, "-gravity", "west", "-append", "-extent", f"{width}x{len(lines) * int(ps * 1.55) + 6}", out)


def info(tag):
    return json.loads((W / f"dp_{tag}.json").read_text())


def fmt_rounds(j):
    names = ["E_P", "E_λ", "E_λ"]
    return " | ".join(f"{n} {r['e_start']:.3g}→{r['e_end']:.3g}" for n, r in zip(names, j["rounds"]))


def caption(j, name):
    r = j["result"]
    return [
        f"{name}: {j['faces_in']:,} → {j['faces_out']:,} faces, {len(j['rounds'])} rounds × {j['rounds'][0]['iterations']} iterations, {float(r['seconds']):.0f} s",
        fmt_rounds(j),
        f"Hausdorff {j['d_haus']:.3f}, volume {j['vol_ratio']:.2f}, folds {j['folds']}, {j['seam_long_components']} seam curves",
    ]


def hero():
    tmp = W / "dev2/hero"
    tmp.mkdir(exist_ok=True)
    cols = []
    S = 760
    for t, name in (("spot_iso", "spot"), ("bob_iso", "bob"), ("blub_iso", "blub")):
        j = info(t)
        cell(P / f"{t}_in.png", tmp / f"{t}_in.png", S, f"uniform input, {j['faces_in']:,} faces")
        cell(P / f"{t}_out.png", tmp / f"{t}_out.png", S, f"result, {j['faces_out']:,} faces, seams in red")
        cell(P / f"{t}_def.png", tmp / f"{t}_def.png", S, "angle defect on the result (blue −, red +, ±0.05 rad)")
        text_block(caption(j, name), S, tmp / f"{t}_txt.png", 15)
        col = tmp / f"{t}_col.png"
        cv(tmp / f"{t}_in.png", tmp / f"{t}_out.png", tmp / f"{t}_def.png", tmp / f"{t}_txt.png",
           "-background", BG, "-append", col)
        cols.append(col)
    cv(*cols, "-background", BG, "+append", "-bordercolor", BG, "-border", "12", M / "developable_hero.png")


def inputs():
    tmp = W / "dev2/inputs"
    tmp.mkdir(exist_ok=True)
    S = 900
    cols = []
    for t, title in (("spot_cage", "subdivision cage as input"), ("spot_iso", "uniform input")):
        j = info(t)
        cell(P / f"{t}_in.png", tmp / f"{t}_in.png", S, f"{title}: {j['faces_in']:,} faces")
        cell(P / f"{t}_out.png", tmp / f"{t}_out.png", S, f"result, {j['faces_out']:,} faces")
        r = j["result"]
        text_block([f"{j['faces_in']:,} → {j['faces_out']:,} faces, same recipe (E_P, then E_λ ×2, 500 iterations each)",
                    fmt_rounds(j),
                    f"Hausdorff to the dense spot {j['d_haus']:.3f} of the diagonal, volume {j['vol_ratio']:.2f}, folds {j['folds']}"],
                   S, tmp / f"{t}_txt.png", 15, bold_first=False)
        col = tmp / f"{t}_col.png"
        cv(tmp / f"{t}_in.png", tmp / f"{t}_out.png", tmp / f"{t}_txt.png", "-background", BG, "-append", col)
        cols.append(col)
    cv(*cols, "-background", BG, "+append", "-bordercolor", BG, "-border", "12", M / "developable_inputs.png")


ENERGIES = [("covariance", "E_λ (covariance)"), ("combinatorial", "E_P (combinatorial)"),
            ("combinatorial-max", "combinatorial, max"), ("covariance-max", "covariance, max"),
            ("intrinsic", "intrinsic"), ("intrinsic-max", "intrinsic, max")]


def energies():
    tmp = W / "dev2/energies"
    tmp.mkdir(exist_ok=True)
    S = 640
    cells = []
    for e, title in ENERGIES:
        j = info(f"spot_e_{e}")
        r = j["rounds"]
        cell(P / f"e_{e}.png", tmp / f"{e}.png", S, title)
        text_block([f"E {r[0]['e_start']:.3g} → {r[-1]['e_end']:.3g}, {j['faces_out']:,} faces, {float(j['result']['seconds']):.0f} s",
                    f"Hausdorff {j['d_haus']:.3f}, volume {j['vol_ratio']:.2f}, folds {j['folds']}, {j['seam_long_components']} seam curves"],
                   S, tmp / f"{e}_txt.png", 14, bold_first=False)
        cv(tmp / f"{e}.png", tmp / f"{e}_txt.png", "-background", BG, "-append", tmp / f"{e}_c.png")
        cells.append(tmp / f"{e}_c.png")
    cv(*cells[:3], "-background", BG, "+append", tmp / "r1.png")
    cv(*cells[3:], "-background", BG, "+append", tmp / "r2.png")
    cv(tmp / "r1.png", tmp / "r2.png", "-background", BG, "-append", "-bordercolor", BG, "-border", "12",
       M / "developable_energies.png")


def crumpled():
    tmp = W / "dev2/crumpled"
    tmp.mkdir(exist_ok=True)
    S = 800
    j = info("crumpled")
    cell(P / "crumpled_in.png", tmp / "in.png", S, f"crumpleddevelopable, {j['faces_in']:,} faces")
    cell(P / "crumpled_out.png", tmp / "out.png", S, f"result, {j['faces_out']:,} faces, seams in red")
    cell(P / "crumpled_def.png", tmp / "def.png", S, "angle defect (blue −, red +, ±0.05 rad)")
    text_block([f"{j['faces_in']:,} → {j['faces_out']:,} faces, the paper's recipe, 3 × 500 iterations, {float(j['result']['seconds']):.0f} s",
                fmt_rounds(j),
                f"edges above 30°: {j['e30']} of {j['edges']}, folds {j['folds']}, {j['seam_long_components']} seam curves, angle defect sum {j['sum_abs_defect']:.2f} rad"],
               3 * S, tmp / "txt.png", 15, bold_first=False)
    cv(tmp / "in.png", tmp / "out.png", tmp / "def.png", "+append", tmp / "row.png")
    cv(tmp / "row.png", tmp / "txt.png", "-background", BG, "-append", "-bordercolor", BG, "-border", "12",
       M / "developable_crumpled.png")


if __name__ == "__main__":
    for k in (sys.argv[1:] or ["all"]):
        for f in ((hero, inputs, energies, crumpled) if k == "all" else (globals()[k],)):
            f()
