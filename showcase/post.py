"""Post-processing: composite renders over the background, add labels, encode videos.

post.py still <in.png> <out.png> [label]
post.py video <frames-dir> <out-basename> [labels.json] [fps] [gif-width]
"""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BG = "#f3f4f6"
FONT = "DejaVu-Sans"
INK = "#1f2329"


def composite(src, dst, label=None, pointsize=None):
    w = int(subprocess.check_output(["identify", "-format", "%w", str(src)]))
    ps = pointsize or max(14, w // 48)
    cmd = ["convert", str(src), "-background", BG, "-alpha", "remove", "-alpha", "off"]
    if label:
        cmd += ["-font", FONT, "-pointsize", str(ps), "-fill", INK, "-gravity", "southwest",
                "-annotate", f"+{ps}+{int(ps * 0.9)}", label]
    cmd.append(str(dst))
    subprocess.run(cmd, check=True)


def video(frames, out, labels=None, fps=30, gif_width=720):
    frames = Path(frames)
    comp = frames.parent / (frames.name + "_comp")
    comp.mkdir(exist_ok=True)
    files = sorted(frames.glob("[0-9]*.png"))
    labs = json.load(open(labels)) if labels else [None] * len(files)
    with ThreadPoolExecutor(32) as ex:
        list(ex.map(lambda i: composite(files[i], comp / files[i].name, labs[i] if i < len(labs) else None),
                    range(len(files))))
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pat = str(comp / "%04d.png")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pat,
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-preset", "slow",
                    "-movflags", "+faststart", str(out.with_suffix(".mp4"))], check=True)
    gfps = min(fps, 20)
    filt = (f"fps={gfps},scale={gif_width}:-1:flags=lanczos,split[a][b];"
            f"[a]palettegen=max_colors=192:stats_mode=diff[p];[b][p]paletteuse=dither=sierra2_4a")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", pat,
                    "-filter_complex", filt, str(out.with_suffix(".gif"))], check=True)
    # Poster frame.
    print("wrote", out.with_suffix(".mp4"), out.with_suffix(".gif"))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "still":
        composite(a[1], a[2], a[3] if len(a) > 3 else None)
    elif a[0] == "video":
        video(a[1], a[2], a[3] if len(a) > 3 and a[3] != "-" else None,
              int(a[4]) if len(a) > 4 else 30, int(a[5]) if len(a) > 5 else 720)
