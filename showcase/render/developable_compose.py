"""Compose the developability stills of media/ from the renders of developable_stills.sh
(work/dev/*.png) and the numbers of the runs (work/developable_*.json). Run:
    uv run --no-project --with numpy python render/developable_compose.py
"""
import json
import subprocess
from pathlib import Path

import numpy as np

SHOW = Path(__file__).resolve().parents[1]
W, T, M = SHOW / "work", SHOW / "work/dev", SHOW / "media"
BG, INK, F = "#f3f4f6", "#1f2329", "DejaVu-Sans"
M.mkdir(exist_ok=True)


def sh(*cmd):
    subprocess.run([str(c) for c in cmd], check=True)


def flat(name):
    out = T / f"flat_{name}.png"
    sh("convert", T / f"{name}.png", "-background", BG, "-alpha", "remove", "-alpha", "off", out)
    return out


def trim_h(img, out, pad=24):
    sh("convert", img, "-fuzz", "1.5%", "-trim", "+repage", "-bordercolor", BG, "-border", pad, out)
    return out


def caption(img, out, title, lines=(), size=34):
    """A title strip above and note lines below the image."""
    title = title.replace("%", "%%")
    lines = [x.replace("%", "%%") for x in lines]
    sh("convert", img, "+repage", "-background", BG, "-gravity", "north", "-splice", f"0x{size * 2}",
       "-font", F, "-pointsize", size, "-fill", INK, "-gravity", "northwest",
       "-annotate", f"+28+{size // 2}", title, out)
    if lines:
        sh("convert", out, "+repage", "-background", BG, "-gravity", "south",
           "-splice", f"0x{int(size * 1.35 * len(lines)) + size // 2}",
           "-font", F, "-pointsize", int(size * 0.8), "-fill", INK, "-gravity", "southwest",
           "-annotate", f"+28+{size // 3}", "\n".join(lines), out)
    return out


def colorbar(path, stops, w=700, h=26, ticks=(), labels_size=24, title=""):
    """A horizontal color bar from hex stops with tick labels [(pos 0..1, text)]."""
    sys_stops = np.array([[int(s[i:i + 2], 16) for i in (1, 3, 5)] for s in stops], dtype=float)
    x = np.linspace(0, 1, len(sys_stops))
    t = np.linspace(0, 1, w)
    rgb = np.stack([np.interp(t, x, sys_stops[:, k]) for k in range(3)], axis=-1).astype(np.uint8)
    ppm = path.with_suffix(".ppm")
    with open(ppm, "wb") as f:
        f.write(f"P6 {w} {h} 255\n".encode())
        f.write(np.repeat(rgb[None], h, axis=0).tobytes())
    cmd = ["convert", ppm, "-bordercolor", BG, "-border", "80x4", "-background", BG,
           "-gravity", "south", "-splice", f"0x{labels_size + 14}", "-font", F,
           "-pointsize", labels_size, "-fill", INK, "-gravity", "southwest"]
    for pos, text in ticks:
        cmd += ["-annotate", f"+{int(80 + pos * w - len(text) * labels_size * 0.3)}+2", text]
    if title:
        cmd += ["-gravity", "north", "-splice", f"0x{labels_size + 10}", "-gravity", "northwest",
                "-annotate", f"+80+4", title]
    cmd.append(path)
    sh(*cmd)
    return path


def stack_bar(img, bar, out):
    sh("convert", img, "+repage", bar, "-background", BG, "-gravity", "west", "-append", out)
    return out


def with_bar(img, bar, out, title, lines=(), size=34):
    """Caption first, then the color bar underneath the notes."""
    caption(img, T / "wb.png", title, lines, size)
    stack_bar(T / "wb.png", bar, out)
    return out


def info(tag):
    return json.loads((W / f"developable_{tag}.json").read_text())


def g(x):
    return f"{x:.3g}"


ENERGY_STOPS = ["#f7f7f4", "#fff2b0", "#fdc761", "#f4823c", "#d03a2f", "#8a1a4a", "#3b0f4f"]
RDBU = ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#e6e4e2", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"]

# --------------------------------------------------------------------------- hero, spot
i = info("f_spot_i6000")
ebar = colorbar(T / "bar_energy.png", ENERGY_STOPS, ticks=[(0, "1e-5"), (0.25, "1e-4"), (0.5, "1e-3"), (0.75, "1e-2"), (1, "1e-1")],
                title="vertex energy E_lambda (log scale, fixed)")
note = [f"spot, {i['n_vertices_in']:,} vertices: {i['iterations']:,} L-BFGS iterations of E_lambda, no subdivision, "
        f"{i['n_vertices_out']:,} vertices after remeshing",
        f"energy {g(i['energy_initial'])} -> {g(i['energy_final'])}, {i['seconds']:.0f} s (one thread); "
        f"edges at |dihedral angle| > 30 degrees: {i['edges_above_30deg_out']} of {i['n_edges_out']}",
        "edge shade: white = flat, dark = dihedral angle of 30 degrees or more"]
caption(trim_h(flat("hero"), T / "t.png"), M / "developable_hero.png",
        "Make developable: input (left) and result (right)", note)
stack_bar(trim_h(flat("energy"), T / "t.png"), ebar, T / "e0.png")
caption(T / "e0.png", M / "developable_energy.png", "Per-vertex energy before and after, one color scale",
        [f"total E_lambda {g(i['energy_initial'])} -> {g(i['energy_final'])}; the energy that is left sits on the creases"])

# --------------------------------------------------------------------------- bob
b = info("f_bob_i4000")
note = [f"bob, {b['n_vertices_in']:,} vertices: {b['iterations']:,} iterations, no subdivision; energy {g(b['energy_initial'])} -> {g(b['energy_final'])}, "
        f"{b['seconds']:.0f} s",
        f"edges above 30 degrees: {b['edges_above_30deg_out']} of {b['n_edges_out']} ({b['n_vertices_out']:,} vertices after remeshing)"]
caption(trim_h(flat("hero_bob"), T / "t.png"), M / "developable_bob.png", "Bob made developable", note)
stack_bar(trim_h(flat("energy_bob"), T / "t.png"), ebar, T / "e1.png")
caption(T / "e1.png", M / "developable_bob_energy.png", "Bob: vertex energy before and after",
        [f"E_lambda {g(b['energy_initial'])} -> {g(b['energy_final'])}"])
b2 = info("f_blub_i3000")
caption(trim_h(flat("hero_blub"), T / "t.png"), M / "developable_blub.png", "Blub made developable",
        [f"blub, {b2['n_vertices_in']:,} vertices: {b2['iterations']:,} iterations; energy {g(b2['energy_initial'])} -> {g(b2['energy_final'])}, {b2['seconds']:.0f} s; "
         f"edges above 30 degrees: {b2['edges_above_30deg_out']} of {b2['n_edges_out']}"])

# --------------------------------------------------------------------------- defect and seams
D = np.load(W / "developable_f_spot_i6000.npz")
ci = json.loads((W / "developable_cut_f_spot_i6000.json").read_text())[1]
dbar = colorbar(T / "bar_defect.png", RDBU, ticks=[(0, "-0.1 rad"), (0.5, "0"), (1, "+0.1 rad")],
                title="angle defect per vertex (fixed scale)")
pair = T / "dpair.png"
sh("convert", flat("defect"), flat("defect_back"), "+append", pair)
stack_bar(trim_h(pair, T / "t.png"), dbar, T / "d0.png")
caption(T / "d0.png", M / "developable_defect.png", "Angle defect and the cut: seams run along the creases",
        [f"spot result, seams chosen with threshold {ci['epsilon']:g}: {ci['n_seam_edges']} seam edges (red), {ci['n_terminals']} terminals (yellow)",
         f"sum of |angle defect| {np.abs(D['out_defect']).sum():.0f} rad over {len(D['out_defect']):,} vertices; front and back views"])

# --------------------------------------------------------------------------- energies compared
names = ["covariance", "combinatorial", "combinatorial_max", "covariance_max", "intrinsic", "intrinsic_max"]
labels = {"covariance": "covariance (E_lambda)", "combinatorial": "combinatorial (E_P)", "combinatorial_max": "combinatorial max",
          "covariance_max": "covariance max", "intrinsic": "intrinsic", "intrinsic_max": "intrinsic max"}
panels = []
p0 = T / "p_in.png"
sh("convert", trim_h(flat("cmp_in"), T / "t.png", 4), "-resize", "x760", p0)
caption(p0, p0, "input", [f"{info('cmp_covariance')['n_vertices_in']:,} vertices"], size=30)
panels.append(p0)
for n in names:
    ii = info(f"cmp_{n}")
    p = T / f"p_{n}.png"
    sh("convert", trim_h(flat(f"cmp_{n}"), T / "t.png", 4), "-resize", "x760", p)
    caption(p, p, labels[n], [f"energy {g(ii['energy_initial'])} -> {g(ii['energy_final'])}", f"edges > 30 deg: {ii['edges_above_30deg_out']}   seams {ii['n_seam_edges']}"], size=30)
    panels.append(p)
# pad to equal sizes and tile 4 x 2
tile = T / "tile.png"
sh("montage", *panels, "-tile", "4x2", "-geometry", "700x700+10+10", "-background", BG, tile)
caption(tile, M / "developable_energies.png", "The six energies, same budget",
        [f"spot ({info('cmp_covariance')['n_vertices_in']:,} vertices), {info('cmp_covariance')['iterations']:,} iterations, no subdivision; energy values are not comparable between panels (different measures)"], size=36)

# --------------------------------------------------------------------------- crumpled
c = info("hero_crumpled")
Dc = np.load(W / "developable_hero_crumpled.npz")
big_in = int((Dc["in_edge_dih"] > np.radians(30)).sum())
note = [f"crumpleddevelopable, {c['n_vertices_in']} vertices, intrinsically flat; {c['subdivisions']} subdivisions x {c['iterations_per_round']} iterations (E_lambda), {c['n_vertices_out']:,} vertices after",
        f"energy {g(c['energy_initial'])} -> {g(c['energy_final'])}, {c['seconds']:.1f} s; edges above 30 degrees: {big_in} of {len(Dc['in_edge_dih'])} -> {c['edges_above_30deg_out']} of {c['n_edges_out']}",
        f"the flow does not preserve isometry: sum of |angle defect| {c['sum_abs_defect_in']:.3g} -> {c['sum_abs_defect_out']:.3g} rad"]
stack_bar(trim_h(flat("crumpled_energy"), T / "t.png"), ebar, T / "e2.png")
sh("convert", trim_h(flat("crumpled"), T / "t.png"), T / "e2.png", "-background", BG, "-append", T / "c0.png")
caption(T / "c0.png", M / "developable_crumpled.png", "Crumpled sheet: wrinkles go, straight rulings appear", note, size=32)

# --------------------------------------------------------------------------- cut and flatten
ls_bar = colorbar(T / "bar_ls.png", RDBU, ticks=[(0, "-0.3"), (0.5, "0"), (1, "+0.3")],
                  title="log scale factor u of the flattening (red = stretched in the plane... see note)")


def cut_panel(tag, k, three, lay, out, title, extra=""):
    ci = json.loads((W / f"developable_cut_{tag}.json").read_text())[k]
    lay_img = T / "lay_t.png"
    sh("convert", flat(lay), "-resize", "x900", lay_img)
    sh("convert", flat(three), "-resize", "x900", T / "three_t.png")
    sh("convert", T / "three_t.png", lay_img, "+append", T / "cp.png")
    stack_bar(trim_h(T / "cp.png", T / "t.png"), ls_bar, T / "cp2.png")
    caption(T / "cp2.png", out, title,
            [f"threshold {ci['epsilon']:g}: {ci['n_seam_edges']} seam edges, {ci['n_terminals']} terminals, {ci['n_charts']} chart, {ci['n_triangles']:,} triangles",
             f"mean |u| {ci['mean_abs_log_scale']:.3g}, quasi-conformal error mean {ci['qc_mean']:.3f} max {ci['qc_max']:.2f}, {ci['flipped']} flipped triangles, overlap {100 * ci['overlap_fraction_max']:.1f} %" + extra])
    return out


ls_bar = colorbar(T / "bar_ls.png", RDBU, ticks=[(0, "-0.3 (shrunk)"), (0.5, "0"), (1, "+0.3 (stretched)")],
                  title="log scale factor u = 0.5 ln(UV area / area), fixed")
cut_panel("f_bob_i4000", 0, "cut_bob0", "lay_bob0", M / "developable_cut_bob.png",
          "Bob: the cut and the flattening (paper Sec. 4.4)")
cut_panel("hero_crumpled", 0, "cut_cr", "lay_cr", M / "developable_cut_crumpled.png",
          "Crumpled sheet after the flow: a nearly isometric flattening")
# epsilon sweep: three columns
cols = []
for k in range(3):
    ci = json.loads((W / "developable_cut_f_bob_i4000.json").read_text())[k]
    sh("convert", flat(f"cut_bob{k}"), "-resize", "x700", T / "s3.png")
    sh("convert", flat(f"lay_bob{k}"), "-resize", "x700", T / "s2.png")
    sh("convert", T / "s3.png", T / "s2.png", "-background", BG, "-append", T / "s1.png")
    caption(T / "s1.png", T / f"sweep{k}.png", f"threshold {ci['epsilon']:g}",
            [f"{ci['n_seam_edges']} seam edges, {ci['n_terminals']} terminals", f"mean |u| {ci['mean_abs_log_scale']:.3g}, {ci['flipped']} flips, overlap {100 * ci['overlap_fraction_max']:.1f} %", f"quasi-conformal max {ci['qc_max']:.2f}"], size=30)
    cols.append(T / f"sweep{k}.png")
sh("convert", *cols, "+append", T / "sw.png")
stack_bar(T / "sw.png", ls_bar, T / "sw2.png")
caption(T / "sw2.png", M / "developable_cut_sweep.png", "Bob: a lower threshold cuts through more vertices and flattens with less distortion", [], size=34)
print("composed")
