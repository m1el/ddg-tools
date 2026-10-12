#!/usr/bin/env python3
"""Stage the Blender extension for `blender --command extension build`.

    blender_stage_extension.py --source blender/ddg_uv_unwrap --out DIR
                               [--license LICENSE] [--placeholder] WHEEL [WHEEL ...]

Copies the add-on sources to DIR, puts the wheels in DIR/wheels, and writes DIR's
blender_manifest.toml: the source manifest with the `platforms` and `wheels` lines filled in from
the wheel file names. With --placeholder the extension's name and tagline say that the bundled
native module contains the stand-in pipeline, so that the built zip cannot be mistaken for the
real unwrap.

Standard library only.
"""

import argparse
import pathlib
import re
import shutil
import sys
import tomllib

MARKER = "# @GENERATED@"


def blender_platforms(wheel_name):
    """Blender platform ids for a wheel file name, from its platform tag."""
    tag = wheel_name[:-len(".whl")].split("-")[-1]
    if tag == "any":
        return []
    if tag.startswith(("manylinux", "musllinux", "linux")):
        system = "linux"
    elif tag.startswith("win"):
        system = "windows"
    elif tag.startswith("macosx"):
        system = "macos"
    else:
        raise SystemExit("unknown platform tag {!r} in {}".format(tag, wheel_name))
    if tag.endswith("universal2"):
        return [system + "-x64", system + "-arm64"]
    if tag.endswith(("x86_64", "amd64")):
        return [system + "-x64"]
    if tag.endswith(("aarch64", "arm64")):
        return [system + "-arm64"]
    raise SystemExit("unknown machine in platform tag {!r} of {}".format(tag, wheel_name))


def toml_list(items):
    return "[" + ", ".join('"{}"'.format(item) for item in items) + "]"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--source", required=True, type=pathlib.Path)
    parser.add_argument("--out", required=True, type=pathlib.Path)
    parser.add_argument("--license", type=pathlib.Path)
    parser.add_argument("--placeholder", action="store_true")
    parser.add_argument("wheels", nargs="+", type=pathlib.Path)
    args = parser.parse_args()

    manifest_source = (args.source / "blender_manifest.toml").read_text(encoding="utf-8")
    if manifest_source.count(MARKER) != 1:
        raise SystemExit("expected exactly one {!r} line in the source manifest".format(MARKER))
    for wheel in args.wheels:
        if not wheel.is_file() or wheel.suffix != ".whl":
            raise SystemExit("not a wheel file: {}".format(wheel))

    if args.out.exists():
        # Only ever delete a directory that is a previous stage (or empty).
        if any(args.out.iterdir()) and not (args.out / "blender_manifest.toml").is_file():
            raise SystemExit("refusing to replace {}: not a staged extension".format(args.out))
        if args.out.resolve() == args.source.resolve():
            raise SystemExit("--out must differ from --source")
        shutil.rmtree(args.out)
    shutil.copytree(
        args.source, args.out,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "wheels", ".*", "*.zip"),
    )
    (args.out / "wheels").mkdir()
    platforms = []
    wheel_paths = []
    for wheel in sorted(args.wheels, key=lambda w: w.name):
        shutil.copy2(wheel, args.out / "wheels" / wheel.name)
        wheel_paths.append("./wheels/" + wheel.name)
        for platform in blender_platforms(wheel.name):
            if platform not in platforms:
                platforms.append(platform)
    if args.license is not None:
        shutil.copy2(args.license, args.out / "LICENSE")

    generated = []
    if platforms:
        generated.append("platforms = " + toml_list(platforms))
    generated.append("wheels = " + toml_list(wheel_paths))
    manifest = manifest_source.replace(MARKER, "\n".join(generated))
    if args.placeholder:
        manifest, n = re.subn(
            r'(?m)^name = "(.*)"$', r'name = "\1 [placeholder build]"', manifest)
        manifest, m = re.subn(
            r'(?m)^tagline = ".*"$',
            'tagline = "PLACEHOLDER BUILD: planar projection, not the real unwrap"', manifest)
        if (n, m) != (1, 1):
            raise SystemExit("could not mark the manifest as a placeholder build")
    parsed = tomllib.loads(manifest)
    for key in ("id", "version", "name", "wheels"):
        if key not in parsed:
            raise SystemExit("staged manifest lacks {!r}".format(key))
    (args.out / "blender_manifest.toml").write_text(manifest, encoding="utf-8")

    print("staged {} {} ({}) in {}".format(parsed["id"], parsed["version"], parsed["name"], args.out))
    for path in wheel_paths:
        print("  wheel", path)
    print("  platforms", ", ".join(platforms) if platforms else "(any)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
