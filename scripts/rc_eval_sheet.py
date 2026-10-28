# /// script
# requires-python = ">=3.10"
# dependencies = ["pillow"]
# ///
"""One still sheet of the 24 grid stills (phase B), 6 x 4.  uv run scripts/rc_eval_sheet.py runs_tmp/rc-eval/runB/videos OUT.png"""
import sys
from pathlib import Path
from PIL import Image, ImageDraw

d = Path(sys.argv[1])
names = [l.strip().split("/")[-1][:-4] for l in open(Path(__file__).with_name("rc-eval-subset.txt"))]
W, H = 450, 467
sheet = Image.new("RGB", (6 * W, 4 * H), "white")
for i, n in enumerate(names):
    im = Image.open(d / f"{n}.png").convert("RGB").resize((W, H), Image.LANCZOS)
    ImageDraw.Draw(im).text((6, 4), n, fill=(0, 0, 0))
    sheet.paste(im, ((i % 6) * W, (i // 6) * H))
sheet.save(sys.argv[2])
print("wrote", sys.argv[2])
