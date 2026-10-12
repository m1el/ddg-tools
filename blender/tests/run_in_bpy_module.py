# SPDX-License-Identifier: MIT
"""Run the add-on tests in Blender-as-a-Python-module (the `bpy` package the Blender Foundation
publishes on PyPI), which is an official build with the Python version official builds bundle.

    python run_in_bpy_module.py EXTENSION.zip [-k NAME]

`scripts/blender-test-official.sh` calls this through `uv run --with bpy==<version> --with numpy`
with the Blender user directories pointed at a scratch directory. The zip is installed with the
same operator `blender --command extension install-file` uses, then `test_addon.py` runs
unchanged. `ddg_tools` must not be installed in the environment: it has to come out of the zip.
"""

import os
import runpy
import sys

# Before Blender starts: it must not look at the real user directory.
for name in ("BLENDER_USER_RESOURCES", "BLENDER_USER_CONFIG", "BLENDER_USER_SCRIPTS",
             "BLENDER_USER_EXTENSIONS"):
    if not os.environ.get(name):
        raise SystemExit("{} is not set; refusing to run against the real user directory".format(
            name))

import bpy  # noqa: E402

zip_path = os.path.abspath(sys.argv[1])
extra = sys.argv[2:]

try:
    import ddg_tools
except ImportError:
    pass
else:
    raise SystemExit("ddg_tools is importable before the extension is installed ({}); "
                     "the test would not prove anything".format(ddg_tools.__file__))

# The module has no --factory-startup; this is its equivalent.
bpy.ops.wm.read_factory_settings()
print("Blender", bpy.app.version_string, "as a module, Python", sys.version.split()[0])
print("user resources:", bpy.utils.resource_path('USER'))
assert os.path.realpath(bpy.utils.resource_path('USER')).startswith(
    os.path.realpath(os.environ["BLENDER_USER_RESOURCES"]))
result = bpy.ops.extensions.package_install_files(
    filepath=zip_path, repo="user_default", enable_on_install=False)
assert result == {'FINISHED'}, result

sys.argv = [sys.argv[0], "--", "--module", "bl_ext.user_default.ddg_uv_unwrap",
            "--zip", zip_path, "--bpy-module"] + extra
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_addon.py"),
               run_name="__main__")
