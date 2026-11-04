"""Montage of the renders (made with OUR crate's renderer): rows = views, columns = truth, N, E-ref, E-geo.
usage: montage.py <render_dir> <case> <out.png> [views...]"""
import sys
from PIL import Image, ImageDraw, ImageFont
rd, case, out = sys.argv[1:4]
views = sys.argv[4:] or ['view', 'cutaway', 'cutplane']
cols = [('truth', 'truth (fine mesh)'), ('N', '(N) network'), ('Eref', '(E-ref) exact, paper rule'), ('Egeo', '(E-geo) exact, Geometric')]
ims = [[Image.open(f'{rd}/{case}_{m}_{v}.png').convert('RGB') for m, _ in cols] for v in views]
w, h = ims[0][0].size
s = min(1.0, 1024 / (4 * w))
tw, th = int(w * s), int(h * s)
head = 22
M = Image.new('RGB', (4 * tw, head + len(views) * th + 18), 'white')
d = ImageDraw.Draw(M)
try:
    f = ImageFont.truetype('DejaVuSans.ttf', 13)
except OSError:
    f = ImageFont.load_default()
for j, (_, lab) in enumerate(cols):
    d.text((j * tw + 6, 4), lab, fill=(11, 11, 11), font=f)
for i, row in enumerate(ims):
    for j, im in enumerate(row):
        M.paste(im.resize((tw, th), Image.LANCZOS), (j * tw, head + i * th))
d.text((6, head + len(views) * th + 2), f'{case}: sphere-traced zero level sets, rendered with the ddg-points-as-tori CPU renderer (ours); same inputs as the authors\' program', fill=(82, 81, 78), font=f)
M.save(out)
print(out, M.size)
